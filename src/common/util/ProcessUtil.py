"""子进程生命周期管理:独立进程组 + 注册表 + 分级终止。

设计约定(维护者请先读这里):

1. CmdUtil 启动的每个外部命令都必须经 launch_managed():子进程自成
   进程组(start_new_session=True),不会随终端 Ctrl+C 直接收到信号——
   信号由父进程(run.py)统一转发,行为与启动方式无关。
2. 不变量:CmdUtil 的启动函数(_run_cmd / _run_cmds_parallel)以
   try/finally 调 terminate_all_children() 收尾,函数返回或抛出时
   受管子进程连同其进程组内成员必已清场。Unix 不提供"父死子亡"
   语义,孤儿清理只能显式做,这里是唯一实现出口。
3. run.py 的退出路径只调 terminate_all_children() 兜底清扫,
   信号处理器里不写任何 kill 逻辑。
4. 分级终止:SIGINT(优雅,snakemake 会停作业、清临时文件)→ 宽限
   DEFAULT_GRACE_S → SIGTERM → 短宽限 → SIGKILL;最后对发过信号的
   进程组做一次 best-effort SIGKILL 组清扫(组内可能有继承 SIG_IGN
   的漏网成员,如非交互 shell 以 & 起的后台任务);force=True 跳过
   优雅阶段直接 SIGKILL(二次信号、兜底清扫用)。
"""
import logging
import os
import signal
import subprocess
import time
from typing import List, Optional

DEFAULT_GRACE_S = 30
KILL_GRACE_S = 5

# 受管子进程注册表(launch_managed 登记,terminate_all_children 清空)。
_children: List[subprocess.Popen] = []


class SignalExit(BaseException):
    """外部终止信号(SIGTERM 等)转成的异常。

    继承 BaseException 而非 Exception:绝不允许被业务代码的
    except Exception 误吞(例如 run.py 测试模式按工作流容错的循环),
    必须一路走到 __main__ 的统一退出路径。
    """

    def __init__(self, signum: int):
        super().__init__(f"terminated by signal {signum}")
        self.signum = signum

    @property
    def exit_code(self) -> int:
        return 128 + self.signum


def launch_managed(cmd: List[str], cwd: Optional[str] = None,
                   env: Optional[dict] = None) -> subprocess.Popen:
    """启动自成进程组的受管子进程并登记注册表。

    stdout/stderr 均为 PIPE + text,与原 subprocess.run(capture_output=True)
    的可见性一致。
    """
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=cwd,
        env=env,
        start_new_session=True,
    )
    _children.append(p)
    return p


def _signal_children(children: List[subprocess.Popen], sig: int) -> None:
    """按进程组给子进程发信号(只波及子进程自己的组,永不误伤父进程)。"""
    for p in children:
        try:
            os.killpg(os.getpgid(p.pid), sig)
        except (ProcessLookupError, PermissionError):
            try:
                p.send_signal(sig)
            except ProcessLookupError:
                pass


def _wait_children(children: List[subprocess.Popen], timeout: float) -> bool:
    """等待 children 全部退出;True=全部已退出, False=超时。"""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if all(p.poll() is not None for p in children):
            return True
        time.sleep(0.2)
    return all(p.poll() is not None for p in children)


def terminate_all_children(grace_s: float = DEFAULT_GRACE_S,
                           force: bool = False,
                           logger: Optional[logging.Logger] = None) -> None:
    """终止所有受管子进程(连同其进程组),幂等,可重复调用。

    force=False:SIGINT → 宽限 grace_s → SIGTERM → 短宽限 → SIGKILL。
    force=True :直接 SIGKILL 并立即返回(不做等待),用于二次信号和
                退出前兜底清扫。
    无论走哪条路,最后都对发过信号的进程组做一次 SIGKILL 组清扫,
    保证"函数退出时整组无生还"(组内可能有 SIGINT 被继承为 SIG_IGN
    的漏网成员,优雅信号清不掉它们)。
    调用结束后注册表清空:子进程的生命周期所有权就此了结。
    """
    log = logger or logging.getLogger(__name__)
    alive = [p for p in _children if p.poll() is None]
    # 记录涉及的进程组(子进程均自成进程组,绝不等于父进程的组)。
    pgids = set()
    for p in alive:
        try:
            pgids.add(os.getpgid(p.pid))
        except ProcessLookupError:
            pass
    if alive:
        if force:
            _signal_children(alive, signal.SIGKILL)
        else:
            log.info("terminating %d child process(es): SIGINT -> %ds grace",
                     len(alive), grace_s)
            _signal_children(alive, signal.SIGINT)
            if not _wait_children(alive, grace_s):
                alive = [p for p in alive if p.poll() is None]
                if alive:
                    log.warning("%d child process(es) ignored SIGINT after %ds; "
                                "escalating to SIGTERM", len(alive), grace_s)
                    _signal_children(alive, signal.SIGTERM)
                    if not _wait_children(alive, KILL_GRACE_S):
                        alive = [p for p in alive if p.poll() is None]
                        if alive:
                            log.warning("%d child process(es) still alive; "
                                        "escalating to SIGKILL", len(alive))
                            _signal_children(alive, signal.SIGKILL)
                            _wait_children(alive, KILL_GRACE_S)
    # 组清扫:登记子进程退出后,其进程组里可能仍有漏网成员(忽略
    # SIGINT 的后台任务等)。best-effort SIGKILL,组已空则无害。
    for pgid in pgids:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
    _children.clear()
