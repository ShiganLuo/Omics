"""ProcessUtil 单测:进程组隔离、优雅转发、分级升级、force、注册表清空。

运行:python3 test/test_process_util.py(或 pytest 收集)。
覆盖的正是 CmdUtil 依赖的核心不变量:terminate_all_children() 之后
不残留任何受管子进程(含其进程组内孙进程)。
"""
import os
import signal
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.common.util.ProcessUtil import (
    SignalExit,
    _children,
    launch_managed,
    terminate_all_children,
)


def _write_script(tmpdir, name, body):
    path = os.path.join(tmpdir, name)
    with open(path, "w") as fh:
        fh.write("#!/bin/bash\n" + body)
    os.chmod(path, 0o755)
    return path


def _alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False


def test_graceful_sigint_forwarded_to_child_group():
    """SIGINT 转发到子进程组:trap 收到并自行退出,孙进程同组被带走。"""
    with tempfile.TemporaryDirectory() as tmpdir:
        marker = os.path.join(tmpdir, "trap.log")
        gpfile = os.path.join(tmpdir, "grandchild.pid")
        script = _write_script(tmpdir, "child.sh", f"""
trap 'echo got-INT >> {marker}; exit 7' INT
sleep 60 &
echo $! > {gpfile}
wait
""")
        p = launch_managed(["bash", script])
        time.sleep(0.5)
        grandchild = int(open(gpfile).read().strip())
        assert _alive(grandchild)

        terminate_all_children(grace_s=5)

        assert p.poll() is not None, "child must be dead"
        assert p.returncode == 7, f"child should exit via its INT trap, got {p.returncode}"
        assert os.path.exists(marker), "SIGINT must reach the child (forwarding broken)"
        # 组清扫 SIGKILL 已投递;僵尸被 init 收走前 kill(pid,0) 仍成功,允许短暂延迟。
        for _ in range(30):
            if not _alive(grandchild):
                break
            time.sleep(0.1)
        assert not _alive(grandchild), "grandchild in the child group must not survive"
        assert not _children, "registry must be empty after terminate"


def test_escalation_to_sigkill():
    """忽略 SIGINT/SIGTERM 的子进程最终被 SIGKILL 收走。"""
    with tempfile.TemporaryDirectory() as tmpdir:
        script = _write_script(tmpdir, "stubborn.sh", "trap '' INT TERM\nsleep 60\n")
        p = launch_managed(["bash", script])
        time.sleep(0.3)

        terminate_all_children(grace_s=1)

        assert p.poll() is not None, "stubborn child must be dead"
        assert p.returncode == -signal.SIGKILL, f"expected SIGKILL, got {p.returncode}"
        assert not _children


def test_force_skips_grace():
    """force=True 立即 SIGKILL,不做等待(二次信号/兜底清扫的逃生口)。"""
    with tempfile.TemporaryDirectory() as tmpdir:
        script = _write_script(tmpdir, "sleeper.sh", "sleep 60\n")
        p = launch_managed(["bash", script])
        time.sleep(0.3)

        t0 = time.monotonic()
        terminate_all_children(force=True)
        elapsed = time.monotonic() - t0

        assert elapsed < 2, f"force must not wait for grace, took {elapsed:.1f}s"
        # force 是 fire-and-forget:SIGKILL 已投递,给内核短暂收尸时间。
        for _ in range(25):
            if p.poll() is not None:
                break
            time.sleep(0.1)
        assert p.poll() is not None
        assert not _children


def test_terminate_is_idempotent_and_ignores_exited_children():
    """重复调用是无害 no-op;已退出子进程不进终止流程。"""
    with tempfile.TemporaryDirectory() as tmpdir:
        script = _write_script(tmpdir, "quick.sh", "exit 0\n")
        p = launch_managed(["bash", script])
        assert p.wait(timeout=5) == 0

        terminate_all_children(grace_s=1)  # 已退出:不应报错、不应等待
        terminate_all_children(grace_s=1)  # 注册表已空:纯 no-op
        assert not _children


def test_signal_exit_semantics():
    e = SignalExit(signal.SIGTERM)
    assert e.signum == signal.SIGTERM
    assert e.exit_code == 143
    assert not isinstance(e, Exception), "SignalExit must not be catchable as Exception"
    assert isinstance(e, BaseException)


def test_group_isolation_parent_untouched():
    """killpg 只打子进程组:父进程(测试自身)不被误伤。"""
    with tempfile.TemporaryDirectory() as tmpdir:
        script = _write_script(tmpdir, "sleeper.sh", "sleep 60\n")
        launch_managed(["bash", script])
        time.sleep(0.3)
        terminate_all_children(grace_s=5)
    assert _alive(os.getpid()), "parent must survive child-group termination"


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                failures += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if failures else 0)
