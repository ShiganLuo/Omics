from snakemake.logging import logger
import time
import os
import sys

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
meta = config.get("meta", "")

# Import serialization utilities
_bin_dir = os.path.join(ROOT_DIR, "modules", "mimseq", "bin")
if _bin_dir not in sys.path:
    sys.path.insert(0, _bin_dir)
from utils import extract_condition

rule mimseq_prepare_sample_data:
    """Prepare sample_data.tsv for mimseq."""
    input:
        meta = meta,
    output:
        sample_data = outdir + "/sample_data.tsv",
    log:
        logdir + "/mimseq_prepare_sample_data.log"
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()

            sample_condition_dict = extract_condition(input.meta)
            os.makedirs(outdir, exist_ok=True)

            with open(output.sample_data, "w") as f:
                for sample_id in samples:
                    if sample_id not in sample_condition_dict:
                        raise ValueError(f"Sample {sample_id} not found in metadata file {input.meta}.")
                    design = sample_condition_dict[sample_id]
                    f.write(f"{os.path.join(indir, sample_id, f'{sample_id}.single.fq.gz')}\t{design}\n")

            with open(log_path, "a") as f:
                f.write(f"Created sample_data.tsv with {len(samples)} samples\n")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"mimseq_prepare_sample_data failed: {e}\n")
            raise

# Include sub-modules
include: "tRNAtools/tRNAtools.smk"
include: "align/align.smk"
include: "clusters/clusters.smk"
include: "mods/mods.smk"
include: "coverage/coverage.smk"
include: "deseq/deseq.smk"

rule mimseq_all:
    """Run all mimseq modules in sequence."""
    input:
        deseq_done = outdir + "/deseq.done",

rule mimseq_result:
    """Aggregate mimseq outputs as dependency endpoint."""
    input:
        outdir = outdir + "/deseq.done",
