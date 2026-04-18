import logging
logger = logging.getLogger(__name__)
aligner = config.get('Procedure',{}).get('aligner')
SOAPnuke_cofig = {
        "indir": workdir,
        "outdir":  f"{outdir}/SOAPnuke",
        "logdir": logdir
}
module SOAPnuke:
    snakefile: "../modules/SOAPnuke/SOAPnuke.smk"
    config: SOAPnuke_cofig
logger.info(f"SOAPnuke config: {SOAPnuke_cofig}")
use rule soapnuke_filter_paired from SOAPnuke as CoCulture_soapnuke_filter_paired
use rule soapnuke_filter_single from SOAPnuke as CoCulture_soapnuke_filter_single

    
hisat2_config = {
        "indir": SOAPnuke_cofig["outdir"],
        "outdir":  f"{outdir}/hisat2",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "hisat2": config.get('Procedure',{}).get('hisat2')
        },
        "genome": {
            "GRCm39": {
                "fasta": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa",
                "index_prefix": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/hisat2/genome"
            },
            "GRCh38": {
                "fasta": "/data/pub/zhousha/Reference/human/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa",
                "index_prefix": None
            }
        }
    }

module hisat2:
    snakefile: "../modules/hisat2/polygenomes/hisat2.smk"
    config: hisat2_config
logger.info(f"hisat2_config: {hisat2_config}")

use rule hisat2_align from hisat2 as CoCulture_hisat2_align
use rule hisat2_index from hisat2 as CoCulture_hisat2_index
