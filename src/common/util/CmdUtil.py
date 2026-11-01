from typing import List, Optional
import subprocess
import shutil
import logging
import shlex
import os
logger = logging.getLogger(__name__)
def _shell_join(cmd: List[str]) -> str:
    """Render a command list as a shell-like string for logging.

    Parameters
    ----------
    cmd : List[str]
        Command tokens.

    Returns
    -------
    str
        A shell-escaped command string compatible with Python 3.6.
    """
    return " ".join(shlex.quote(part) for part in cmd)
def _run_cmd_p36(cmd: List[str]) -> str:
    """
    Execute complex command and return stdout. python3.6

    - Command not found: give clear message
    - Command execution failed: print stdout/stderr

    Parameters
    ----------
    cmd : List[str]
        Command and arguments, e.g., ["ls", "-l"]

    Returns
    -------
    str
        Standard output of the command.

    Raises
    ------
    RuntimeError
        If command not found or execution fails.
    """
    cmd_str = _shell_join(cmd)
    cmd_bin = cmd[0]

    logger.info(f"Running: {cmd_str}")

    # Precheck: is command available?
    if shutil.which(cmd_bin) is None:
        logger.error(f"Command not found: '{cmd_bin}'")
        logger.error("Please make sure it is installed and in $PATH")
        raise RuntimeError(f"Command not found: {cmd_bin}")

    try:
        result = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout_bytes, stderr_bytes = result.communicate()
        retcode = result.returncode

        # decode bytes to str
        stdout = stdout_bytes.decode("utf-8") if stdout_bytes else ""
        stderr = stderr_bytes.decode("utf-8") if stderr_bytes else ""

        if retcode != 0:
            logger.error(f"Command failed with return code {retcode}")
            logger.error(f"STDOUT:\n{stdout or '[empty]'}")
            logger.error(f"STDERR:\n{stderr or '[empty]'}")
            raise RuntimeError(f"Command execution failed: {cmd_str}")

        if stdout:
            logger.info(f"Command Output:\n{stdout}")

        return stdout

    except OSError as e:
        logger.error(f"Execution failed: {str(e)}")
        raise RuntimeError(f"Command execution failed: {cmd_str}") from e

def _sanitized_env(cmd: List[str]) -> dict:
    """Build a subprocess environment with TMPDIR forced to the workflow's
    own tmp_dir.

    The calling shell may export a host-only TMPDIR (agent scratch dirs etc.)
    that is NOT bind-mounted into apptainer containers — tools inheriting it
    (GNU parallel, fibertools temp files, ...) die at startup inside the
    container. The workflow's own tmp_dir is always bind-mounted (and mapped
    to /tmp inside the container), so force that instead. Relative tmp_dir
    values resolve against the configfile's directory, matching
    run.build_snakemake_cmd.
    """
    env = dict(os.environ)
    cfg_tmp = None
    if "--configfile" in cmd:
        import json as _json
        idx = cmd.index("--configfile")
        if idx + 1 < len(cmd):
            cfg_path = cmd[idx + 1]
            try:
                with open(cfg_path, "r", encoding="utf-8") as _f:
                    _cfg = _json.load(_f)
                cfg_tmp = _cfg.get("tmp_dir")
            except Exception:
                cfg_tmp = None
            if cfg_tmp and not os.path.isabs(cfg_tmp):
                cfg_tmp = os.path.join(os.path.dirname(cfg_path), cfg_tmp)
    if cfg_tmp:
        os.makedirs(cfg_tmp, exist_ok=True)
        env["TMPDIR"] = cfg_tmp
    return env

def _run_cmd(cmd:List, cwd: Optional[str] = None) -> str:
    """
    execute complex command and return stdout
    - Command not found: give clear message
    - Command execution failed: print stdout/stderr
    - cwd: working directory for the subprocess (default: inherit from parent)
    """
    cmd_str = " ".join(cmd)
    cmd_bin = cmd[0]

    logger.info(f"Running: {cmd_str}")

    # precheck：is command available?
    if shutil.which(cmd_bin) is None:
        logger.error(f"Command not found: '{cmd_bin}'")
        logger.error("Please make sure it is installed and in $PATH")
        raise RuntimeError(f"Command not found: {cmd_bin}")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            cwd=cwd,
            env=_sanitized_env(cmd)
        )

        if result.stdout:
            logger.info(f"Command Output:\n{result.stdout}")

        return result.stdout

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with return code {e.returncode}")
        logger.error(f"STDOUT:\n{e.stdout or '[empty]'}")
        logger.error(f"STDERR:\n{e.stderr or '[empty]'}")
        raise RuntimeError(
            f"Command execution failed: {cmd_str}"
        ) from e


def _run_cmds_parallel(cmds: List) -> List[str]:
    """Execute multiple commands in parallel, each with its own cwd.

    Parameters
    ----------
    cmds : list of (cmd_list, cwd) tuples
        Each entry is a (command_tokens, working_directory) pair.

    Returns
    -------
    list of str
        Stdout from each command, in the same order as input.

    Raises
    ------
    RuntimeError
        If any command fails (after all have finished).
    """
    processes: List[subprocess.Popen] = []
    for cmd, cwd in cmds:
        cmd_bin = cmd[0]
        if shutil.which(cmd_bin) is None:
            raise RuntimeError(f"Command not found: {cmd_bin}")
        logger.info(f"Starting (cwd={cwd}): {' '.join(cmd)}")
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=_sanitized_env(cmd),
        )
        processes.append(p)

    results: List[str] = []
    failed: List[tuple[int, str, str]] = []
    for i, p in enumerate(processes):
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            failed.append((i, stdout or "", stderr or ""))
        else:
            if stdout:
                logger.info(f"Process {i} output:\n{stdout}")
            results.append(stdout)

    if failed:
        for idx, stdout, stderr in failed:
            logger.error(f"Process {idx} failed (rc={processes[idx].returncode})")
            logger.error(f"STDOUT:\n{stdout or '[empty]'}")
            logger.error(f"STDERR:\n{stderr or '[empty]'}")
        raise RuntimeError(
            f"{len(failed)}/{len(processes)} parallel commands failed"
        )

    return results