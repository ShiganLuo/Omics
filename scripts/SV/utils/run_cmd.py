import subprocess
import shutil
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)

def run_cmd_list(cmd:list):
        """
        执行外部命令，返回 stdout
        - 命令不存在：给出清晰提示
        - 命令执行失败：打印 stdout / stderr
        """
        cmd_str = " ".join(cmd)
        cmd_bin = cmd[0]

        logger.info(f"Running: {cmd_str}")

        # 1️⃣ 预检查：命令是否存在（比 FileNotFoundError 更友好）
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
