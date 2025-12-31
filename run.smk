shell.prefix("set -x; set -e;")
import logging
import os
from itertools import chain
import sys
from snakemake.io import glob_wildcards
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    stream=sys.stdout,
	datefmt='%Y-%m-%d %H:%M:%S'
)
# containerize: "quay.nju.edu.cn"
EXECUTION_DIR = os.getcwd()
SNAKEFILE_FULL_PATH = workflow.snakefile
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE_FULL_PATH)
indir = config.get('indir', 'data')
outdir = config.get('outdir', 'output')
metadata = config.get('metadata')
logging.info("Workflow RNA-SNP started.")
logging.info(f"metadata file path: {metadata}")
logging.info(f"Input directory: {indir}")
logging.info(f"Output directory: {outdir}")
logging.info(f"Snakefile path: {SNAKEFILE_FULL_PATH}")
logging.info(f"Execution directory: {EXECUTION_DIR}")

sys.path.append(f"{SNAKEFILE_DIR}/utils")
from fastq_utils import MetadataUtils
from smk_utils import get_yaml_path
from typing import DefaultDict, List, Dict
configfilePath = os.path.join(SNAKEFILE_DIR,"config","run.yaml")
configfile: configfilePath
logging.info(f"add cofigfile {configfilePath}")
metadataUtils = MetadataUtils(metadata,indir,f"{outdir}/log/utils/fastq_utils.log")
groups = metadataUtils.run()

# control variable
outfiles = []
genomes = []
paired_samples = []
single_samples = []
all_samples = []
XenofilterR_target_samples = []
single_sample_genome_pairs = [] # [(SE,genome)……]
paired_sample_genome_pairs = [] # [(PE,genome)……]
XenofilterR_target_genome = "human"
XenofilterR_pollution_source_genome = "mouse"
# get outfiles and fill control variable
def get_output_Count(groups:Dict[str, Dict[str, List[str]]]):
    include: "subworkflow/Align/Align.smk"
    include: "subworkflow/Count/Count.smk"

    for genome, library_sample in groups.items():
        genomes.append(genome)
        for libraryStrategy, samples in library_sample.items():
            if libraryStrategy == "PAIRED":
                outfiles.append(f"{outdir}/counts/featureCounts/{genome}/{genome}_paired_count.tsv")
                for sample_id in samples:
                    paired_samples.append(sample_id)
                    all_samples.append(sample_id)
                    paired_sample_genome_pairs.append((sample_id,genome))
                    outfiles.append(f"{outdir}/Align/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
            elif libraryStrategy == "SINGLE":
                outfiles.append(f"{outdir}/counts/featureCounts/{genome}/{genome}_single_count.tsv")
                for sample_id in samples:
                    paired_samples.append(sample_id)
                    all_samples.append(sample_id)
                    single_sample_genome_pairs.append((sample_id,genome))
                    outfiles.append(f"{outdir}/Align/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
            else:
                continue
# get_output_Count(groups)

def get_output_TEtranscripts(groups:Dict[str, Dict[str, List[str]]]):
    include: "subworkflow/Align/Align.smk"
    include: "subworkflow/TEtranscripts/TEtranscripts.smk"
    
    for genome, library_sample in groups.items():
        genomes.append(genome)
        outfiles.append(f"{outdir}/TEtranscripts/TEcount/{genome}/all_TEcount.cntTable")
        for libraryStrategy, samples in library_sample.items():
            if libraryStrategy == "PAIRED":
                for sample_id in samples:
                    paired_samples.append(sample_id)
                    all_samples.append(sample_id)
                    paired_sample_genome_pairs.append((sample_id,genome))
                    outfiles.append(f"{outdir}/TEtranscripts/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
            elif libraryStrategy == "SINGLE":
                outfiles.append(f"{outdir}/counts/featureCounts/{genome}/{genome}_single_count.tsv")
                for sample_id in samples:
                    paired_samples.append(sample_id)
                    all_samples.append(sample_id)
                    single_sample_genome_pairs.append((sample_id,genome))
                    outfiles.append(f"{outdir}/TEtranscripts/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
            else:
                continue    
get_output_TEtranscripts(groups)

def get_output_SNP(groups:Dict[str, Dict[str, List[str]]]):
    include: "subworkflow/SNP/SNP.smk"
    logging.info(f"Include SNP workflow: {SNPSmk}")
def get_output_XenofilterR(groups:Dict[str, Dict[str, List[str]]]):
    include: "subworkflow/XenofilterR/XenofilterR.smk"

logging.info(f"genomes:{genomes}\npaired_sampes:{paired_samples}\nsingle_samples:{single_samples}\nall output files:{outfiles}\n\
XenofilterR_target_samples:{XenofilterR_target_samples}\nXenofilterR_pollution_source_genome:{XenofilterR_pollution_source_genome}\n\
paired_sample_genome_pairs:{paired_sample_genome_pairs}\nsingle_sample_genome_pairs:{single_sample_genome_pairs}")

rule all:
    input:
        outfiles
        # outdir + "/multiqc/multiqc_report.html",


# rule multiqc:
#     input: 
#         get_multiqc_file(paired_samples,single_samples,all_samples,genomes)
#     output:
#         outdir + "/multiqc/multiqc_report.html"
#     params:
#         multiqc_indir = outdir,
#         multiqc_outdir = outdir + "/multiqc/"
#     conda:
#         config['conda']['run']
#     shell:
#         """
#         multiqc {params.multiqc_indir} -o {params.multiqc_outdir}
#         """
    




