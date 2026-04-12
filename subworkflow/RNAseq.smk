import logging
logger = logging.getLogger(__name__)
aligner = config.get('aligner', 'hisat2')
cutadapt_config = {
        "indir": workdir,
        "outdir":  f"{outdir}/cutadapt",
        "logdir": logdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        }
    }
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
logger.info(f"cutadapt_config: {cutadapt_config}")
use rule trimming_Paired from cutadapt as MERIP_trimming_Paired
use rule trimming_Single from cutadapt as MERIP_trimming_Single

if aligner == 'hisat2':
    hisat2_config = {
            "indir": cutadapt_config["outdir"],
            "outdir":  f"{outdir}/hisat2",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "hisat2": config.get('Procedure',{}).get('hisat2')
            },
            "genome": {
                "fasta": config.get('genome',{}).get('fasta')
            }
        }
    module hisat2:
        snakefile: "../modules/hisat2/hisat2.smk"
        config: hisat2_config
    logger.info(f"hisat2_config: {hisat2_config}")
    use rule hisat2_align from hisat2 as MERIP_hisat2_align
    use rule hisat2_index from hisat2 as MERIP_hisat2_index
elif aligner == 'star':
    star_config = {
            "indir": cutadapt_config["outdir"],
            "outdir":  f"{outdir}/star",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "star": config.get('Procedure',{}).get('star')
            },
            "genome": {
                "fasta": config.get('genome',{}).get('fasta')
            }
        }
    module star:
        snakefile: "../modules/star/star.smk"
        config: star_config
    logger.info(f"star_config: {star_config}")
    use rule star_align from star as MERIP_star_align
    use rule star_index from star as MERIP_star_index
else:
    raise ValueError(f"Unsupported aligner: {aligner}")


