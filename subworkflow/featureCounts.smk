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
        snakefile: "../modules/hisat2/hisat2.smk"
        config: hisat2_config
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
        snakefile: "../modules/star/star.smk"
        config: star_config
    use rule star_align from star as RNA_SNP_star_align
else:
    raise ValueError(f"Unsupported aligner: {aligner}. Please choose 'hisat2' or 'star'.")

featureCounts_config = {
        "indir": cutadapt_config["outdir"],
        "outdir":  outdir,
        "Procedure": {
            "featureCounts": config.get('Procedure',{}).get('featureCounts')
        }
    }
module featureCounts:
    snakefile: "../modules/featureCounts/featureCounts.smk"
    config: featureCounts_config
use rule featureCounts from featureCounts as RNA_SNP_featureCounts

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
                    outfiles.append(f"{outdir}/Align/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam")
            elif libraryStrategy == "SINGLE":
                outfiles.append(f"{outdir}/counts/featureCounts/{genome}/{genome}_single_count.tsv")
                for sample_id in samples:
                    paired_samples.append(sample_id)
                    all_samples.append(sample_id)
                    single_sample_genome_pairs.append((sample_id,genome))
                    outfiles.append(f"{outdir}/Align/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam")
            else:
                continue