shell.prefix("set -x; set -e;")
import logging
import os
from snakemake.io import glob_wildcards
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	datefmt='%Y-%m-%d %H:%M:%S'
)
# containerize: "quay.nju.edu.cn"
EXECUTION_DIR = os.getcwd()
SNAKEFILE_FULL_PATH = workflow.snakefile
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE_FULL_PATH)
indir = config.get('indir', '../data')
outdir = config.get('outdir', '../output')
genomes = config.get('genomes', ['human'])
logging.info("Workflow RNA-SNP started.")
logging.info(f"Genomes to be processed: {genomes}")
logging.info(f"Input directory: {indir}")
logging.info(f"Output directory: {outdir}")
logging.info(f"Snakefile path: {SNAKEFILE_FULL_PATH}")
logging.info(f"Execution directory: {EXECUTION_DIR}")

def get_samples(indir:str,outdir:str)->list:
    paired_samples = glob_wildcards(indir + "/{sample_id}_1.fastq.gz").sample_id
    Allsamples = glob_wildcards(indir + "/{sample_id}.fastq.gz").sample_id
    single_samples = [sample for sample in Allsamples if re.match(r"^SRR\d+$", sample)]
    samples = paired_samples + single_samples
    logging.info(f"Detected paired samples: {paired_samples}")
    logging.info(f"Detected single samples: {single_samples}")
    return paired_samples, single_samples, samples
logging.info("preparing samples...")
paired_samples, single_samples, samples = get_samples(indir,outdir)
logging.info("samples prepared.")

def get_snakefile_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the core_snakefile_path/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR, "subworkflow",module_name, f"{module_name}.smk")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module snakefile {module_name}.smk not found at {module_path}")
    return module_path
alignSmk = get_snakefile_path("Align")
include: alignSmk
logging.info(f"Include Align workflow: {alignSmk}")
samples = samples[0:1]  # For test purpose, process only first 1 samples
rule all:
    input:
        expand(outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam", sample_id=samples,genome=genomes)







