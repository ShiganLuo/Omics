import logging
logger = logging.getLogger(__name__)
aligner = config.get('Procedure',{}).get('aligner')
cutadapt_config = {
        "indir": workdir,
        "outdir":  outdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        }
    }
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: config
use rule trimming_Paired from cutadapt as RNA_SNP_trimming_Paired
use rule trimming_Single from cutadapt as RNA_SNP_trimming_Single

if aligner == "hisat2":
    hisat2_config = {
        "indir": cutadapt_config["outdir"],
        "outdir":  outdir,
        "Procedure": {
            "hisat2": config.get('Procedure',{}).get('hisat2')
        }
    }
    module hisat2:
        snakefile: "../modules/hisat2/TEtranscripts/hisat2.smk"
        config: hisat2_config
    use rule hisat2_index from hisat2 as RNA_SNP_hisat2_index
    use rule hisat2_align from hisat2 as RNA_SNP_hisat2_align
elif aligner == "star":
    star_config = {
        "indir": cutadapt_config["outdir"],
        "outdir":  outdir,
        "Procedure": {
            "star": config.get('Procedure',{}).get('star')
        }
    }
    module star:
        snakefile: "../modules/star/TEtranscripts/star.smk"
        config: star_config
    use rule star_align from star as RNA_SNP_star_align
else:
    raise ValueError(f"Unsupported aligner: {aligner}. Please choose 'hisat2' or 'star'.")

def get_output_TEtranscripts(groups:Dict[str, Dict[str, List[str]]]):
    include: "subworkflow/Align/Align.smk"
    include: "subworkflow/TEtranscripts/TEtranscripts.smk"
    
    for genome, library_sample in groups.items():
        genomes.append(genome)
        outfiles.append(f"{outdir}/TEtranscripts/TEcount/{genome}/all_TEcount.tsv")
        for libraryStrategy, samples in library_sample.items():
            if libraryStrategy == "PAIRED":
                for sample_id in samples:
                    paired_samples.append(sample_id)
                    all_samples.append(sample_id)
                    paired_sample_genome_pairs.append((sample_id,genome))
                    outfiles.append(f"{outdir}/TEtranscripts/bam/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam")
            elif libraryStrategy == "SINGLE":
                for sample_id in samples:
                    single_samples.append(sample_id)
                    all_samples.append(sample_id)
                    single_sample_genome_pairs.append((sample_id,genome))
                    outfiles.append(f"{outdir}/TEtranscripts/bam/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam")
            else:
                continue    