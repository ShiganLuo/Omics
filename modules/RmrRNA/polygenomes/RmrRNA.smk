include: "../../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
def get_input_for_extract_rRNA(wildcards):
    """Dynamically determines the input fasta and gtf files for extract_rRNA based on the genome."""
    logger.info(f"[get_input_for_extract_rRNA] called with wildcards: {wildcards}")
    in_dict = {}
    fasta = config.get('genomes', {}).get('references', {}).get(wildcards.genome, {}).get('fasta')
    gtf = config.get('genomes', {}).get('references', {}).get(wildcards.genome, {}).get('gtf')
    if not fasta or not os.path.exists(fasta):
        logger.error(f"Fasta file for genome {wildcards.genome} not found in config or does not exist")
        raise ValueError(f"Fasta file for genome {wildcards.genome} not found in config or does not exist")
    if not gtf or not os.path.exists(gtf):
        logger.error(f"GTF file for genome {wildcards.genome} not found in config or does not exist")
        raise ValueError(f"GTF file for genome {wildcards.genome} not found in config or does not exist")
    in_dict['fasta'] = fasta
    in_dict['gtf'] = gtf
    return in_dict
rule extract_rRNA:
    input:
        unpack(get_input_for_extract_rRNA)
    output:
        rRNA_fasta = outdir + "/{genome}/rRNA.fasta"
    log:
        logdir + "/{genome}/RmrRNA/extract_rRNA.log"
    threads: 2
    params:
        extract_rRNA_script = os.path.join(ROOT_DIR, "modules/RmrRNA/bin/extract_rRNA.py")
    conda:
        "../RmrRNA.yaml"
    container:
        sif("../RmrRNA.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path,'w').close()
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = f"{outdir}/{wildcards.genome}/extract_rRNA.{current_time}.sh"
            cmd = ["python", params.extract_rRNA_script, 
                    "--fasta", input.fasta, 
                    "--gtf", input.gtf, 
                    "--output", output.rRNA_fasta,
                    "--threads", str(threads)
                    ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'extract_rRNA for genome {wildcards.genome} completed'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during extract_rRNA for genome {wildcards.genome}: {e}\n")
            logger.error(f"Error occurred during extract_rRNA for genome {wildcards.genome}: {e}")
            raise e
