import subprocess
import shutil
from pathlib import Path


def gsm_run_info(gsm: str, outfile: str) -> None:
    """
    Query SRA RunInfo for a GSM using NCBI EDirect tools.

    Equivalent to:
        esearch -db gds -query GSM |
        elink -target sra |
        efetch -format runinfo

    :param gsm: GSM ID (e.g. GSM123456)
    :param outfile: Output CSV file
    """
    if not gsm or not outfile:
        raise ValueError("Usage: gsm_run_info(gsm, outfile)")

    print(f"Querying GSM: {gsm}...")

    outfile = Path(outfile)
    tmp_file = outfile.with_suffix(".tmp")

    # 构建 pipeline
    p1 = subprocess.Popen(
        ["esearch", "-db", "gds", "-query", gsm],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    p2 = subprocess.Popen(
        ["elink", "-target", "sra"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    p3 = subprocess.Popen(
        ["efetch", "-format", "runinfo"],
        stdin=p2.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    stdout, stderr = p3.communicate()

    if p3.returncode != 0 or not stdout.strip():
        raise RuntimeError(f"Failed to fetch RunInfo for {gsm}:\n{stderr}")

    tmp_file.write_text(stdout, encoding="utf-8")

    lines = tmp_file.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise RuntimeError(f"No RunInfo returned for {gsm}")

    header = lines[0]
    matched = [line for line in lines[1:] if f",{gsm}," in f",{line},"]

    if not matched:
        print(f"Warning: no RunInfo rows matched GSM {gsm}")

    # outfile 不存在 → 写 header
    if not outfile.exists():
        outfile.parent.mkdir(parents=True, exist_ok=True)
        outfile.write_text(header + "\n", encoding="utf-8")

    with outfile.open("a", encoding="utf-8") as f:
        for line in matched:
            f.write(line + "\n")

    tmp_file.unlink(missing_ok=True)
    print(f"Results saved to {outfile}")

if __name__ == "__main__":
    gsm_run_info("GSM8426614", "data/gsm_runinfo.csv")