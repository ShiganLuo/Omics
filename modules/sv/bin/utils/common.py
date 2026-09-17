from typing import Optional
import logging
logger = logging.getLogger(__name__)
def setup_logger(logger_name: str, level: int = logging.INFO, log_file: Optional[str] = None) -> logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger

def _run_cmd(cmd:List, cwd: Optional[str] = None,logger: logging.Logger = logger) -> str:
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