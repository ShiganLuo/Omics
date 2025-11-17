shell.prefix("set -x; set -e;")
import logging
import os
from itertools import chain
import sys
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
metadata = config.get('metadata')
logging.info("Workflow RNA-SNP started.")
logging.info(f"metadata file path: {metadata}")
logging.info(f"Input directory: {indir}")
logging.info(f"Output directory: {outdir}")
logging.info(f"Snakefile path: {SNAKEFILE_FULL_PATH}")
logging.info(f"Execution directory: {EXECUTION_DIR}")
sys.path.append(f"{SNAKEFILE_DIR}/utils")
from fastq_utils import SNPMetadata
snpMetadata = SNPMetadata(metadata,indir,f"{outdir}/log/utils/fastq_utils.log")
groups = snpMetadata.run()
def get_output_files(groups):
    outfiles = []
    paired_samples = []
    single_samples = []
    all_samples = []
    genomes = []
    for organism, types in groups.items():
        genomes.append(organism)
        for TYPE, samples in types.items():
            if TYPE == "PAIRED":
                paired_samples.append(samples)
                all_samples.append(samples)
            elif TYPE == "SINGLE":
                single_samples.append(samples)
                all_samples.append(samples)
            else:
                continue
    # flatten
    
    paired_samples = list(chain.from_iterable(paired_samples))
    single_samples = list(chain.from_iterable(single_samples))
    all_samples = list(chain.from_iterable(all_samples))
    outfiles = expand(outdir + "/SNP/vcf/filter/{genome}/{sample_id}.vcf.gz",genome=genomes,sample_id=all_samples)
    outfiles_flatten = outfiles
    return outfiles_flatten,paired_samples,single_samples,all_samples,genomes

outfiles_flatten,paired_samples,single_samples,all_samples,genomes = get_output_files(groups)
XenofilterR_target_samples=all_samples
XenofilterR_pollution_source_genome="mouse"
logging.info(f"genomes:{genomes}\npaired_sampes:{paired_samples}\nsingle_samples:{single_samples}\nrule all input files:{outfiles_flatten}")
def get_multiqc_file(paired_samples,single_samples,all_samples,genomes):
    outfiles = []
    outfiles.append(
        expand(
            outdir + "/2pass/{sample_id}/{genome}/{sample_id}Log.final.out",
            sample_id=all_samples,
            genome=genomes
        )
    )
    outfiles.append(
        expand(outdir + "/log/{sample_id}/trimming_statistics_1.txt",sample_id=paired_samples)
    )
    outfiles.append(
        expand(outdir + "/log/{sample_id}/trimming_statistics_2.txt",sample_id=paired_samples)
    )
    outfiles.append(
        expand(
            outdir + "/2pass/{sample_id}/{genome}/{sample_id}Log.final.out",
            sample_id=all_samples,
            genome=genomes
        )
    )
    outfiles.append(
        expand(outdir + "/log/{sample_id}/trimming_statistics.txt",sample_id=single_samples),
    )
    # flatten
    outfiles_flatten = list(chain.from_iterable(outfiles))
    return outfiles_flatten

configfilePath = os.path.join(SNAKEFILE_DIR,"config","run.yaml")
configfile: configfilePath
logging.info(f"add cofigfile {configfilePath}")

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
TEtranscriptsSmk = get_snakefile_path("TEtranscripts")
include: TEtranscriptsSmk
logging.info(f"Include TEtranscripts workflow: {TEtranscriptsSmk}")
SNPSmk = get_snakefile_path("SNP")
include: SNPSmk
logging.info(f"Include SNP workflow: {SNPSmk}")
XenofilterRSmk = get_snakefile_path("XenofilterR")
include: XenofilterRSmk
logging.info(f"Include XenofilterR workflow: {XenofilterRSmk}")


rule all:
    input:
        outfiles_flatten
        # outdir + "/multiqc/multiqc_report.html",


rule multiqc:
    input: 
        get_multiqc_file(paired_samples,single_samples,all_samples,genomes)
    output:
        outdir + "/multiqc/multiqc_report.html"
    params:
        multiqc_indir = outdir,
        multiqc_outdir = outdir + "/multiqc/"
    conda:
        config['conda']['run']
    shell:
        """
        multiqc {params.multiqc_indir} -o {params.multiqc_outdir}
        """
    