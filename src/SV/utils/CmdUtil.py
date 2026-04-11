import subprocess
import shutil
import logging
from typing import List, Optional
import shlex
try:
    from .LogUtil import setup_logger
except Exception:
    from LogUtil import setup_logger

def _run_cmd(cmd:List[str],logger:Optional[logging.Logger] = None) -> str:
    """
    execute complex command and return stdout
    - Command not found: give clear message
    - Command execution failed: print stdout/stderr
    """
    if logger is None:
        logger = setup_logger(__name__)
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
            text=True
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

def _run_cmd_3p6(cmd: List[str], logger: Optional[logging.Logger] = None) -> str:
    """
    Execute complex command and return stdout. python3.6

    - Command not found: give clear message
    - Command execution failed: print stdout/stderr

    Parameters
    ----------
    cmd : List[str]
        Command and arguments, e.g., ["ls", "-l"]
    logger : logging.Logger
        Logger object for logging info/errors

    Returns
    -------
    str
        Standard output of the command.

    Raises
    ------
    RuntimeError
        If command not found or execution fails.
    """
    if logger is None:
        logger = setup_logger(__name__)
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