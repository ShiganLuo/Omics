from typing import List, Optional
import subprocess
import shutil
import logging
import shlex
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
            cwd=cwd
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