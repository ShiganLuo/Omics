import subprocess
import tempfile
import os
import glob
import pandas as pd


def count_terra_reads_pe_seqkit(
    fastq_r1: str,
    fastq_r2: str,
    motif_regex: str = r"(TTAGGG){3,}|(CCCTAA){3,}",
) -> dict:
    """
    Count TERRA-like read pairs from paired-end FASTQ using seqkit.

    A read pair is counted once if either R1 or R2 matches
    telomeric repeats (TTAGGG or CCCTAA).

    Returns
    -------
    dict
        {
            'sample': str,
            'terra_pairs': int,
            'total_pairs': int,
            'terra_ratio': float
        }
    """

    sample = os.path.basename(fastq_r1).replace("_1.fastq.gz", "")

    with tempfile.TemporaryDirectory() as tmp:

        r1_ids = os.path.join(tmp, "r1.ids")
        r2_ids = os.path.join(tmp, "r2.ids")

        # 1️⃣ extract read IDs with telomeric motif
        cmd_r1 = (
            f"seqkit grep -s -r -p '{motif_regex}' {fastq_r1} "
            f"| seqkit seq -n > {r1_ids}"
        )
        cmd_r2 = (
            f"seqkit grep -s -r -p '{motif_regex}' {fastq_r2} "
            f"| seqkit seq -n > {r2_ids}"
        )

        subprocess.run(cmd_r1, shell=True, check=True)
        subprocess.run(cmd_r2, shell=True, check=True)

        # 2️⃣ union read IDs (read-pair level)
        cmd_union = f"cat {r1_ids} {r2_ids} | sort -u | wc -l"
        terra_pairs = int(
            subprocess.check_output(cmd_union, shell=True).strip()
        )

    # 3️⃣ total read pairs (use R1 only)
    cmd_total = ["seqkit", "stats", "-T", fastq_r1]
    out_total = subprocess.check_output(cmd_total, text=True)
    total_pairs = int(out_total.strip().splitlines()[1].split("\t")[3])

    terra_ratio = terra_pairs / total_pairs if total_pairs > 0 else 0.0

    return {
        "sample": sample,
        "terra_pairs": terra_pairs,
        "total_pairs": total_pairs,
        "terra_ratio": terra_ratio,
    }


if __name__ == "__main__":

    results = []

    # 只遍历 R1
    for r1 in sorted(glob.glob("/disk5/luosg/Totipotent20251031/data/fq/*_1.fastq.gz")):
        r2 = r1.replace("_1.fastq.gz", "_2.fastq.gz")

        if not os.path.exists(r2):
            print(f"[WARN] Missing R2 for {r1}, skipped")
            continue

        res = count_terra_reads_pe_seqkit(r1, r2)
        results.append(res)

    df = pd.DataFrame(results)
    df.to_csv("TERRA_from_RNAseq.tsv", sep="\t", index=False)
