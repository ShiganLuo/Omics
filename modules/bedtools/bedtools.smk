include: "../common/common.smk"
from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule iCLIP_bedtools:
    input:
        bam = indir + "/{sample_id}.dedup.bam",
        bai = indir + "/{sample_id}.dedup.bam.bai",
        chromosome_sizes = config.get('genome',{}).get('chrom_sizes')
    output:
        bed = outdir + "/{sample_id}/{sample_id}.bed",
        plus_bedgraph = outdir + "/{sample_id}/{sample_id}.plus.bw",
        minus_bedgraph = outdir + "/{sample_id}/{sample_id}.minus.bw"
    log:
        log = logdir + "/{sample_id}/bedtools.log"
    threads: 4
    conda:
        "bedtools.yaml"
    params:
        bedtools = config.get('Procedure',{}).get('bedtools') or 'bedtools',
        bedGraphToBigWig = config.get('Procedure',{}).get('bedGraphToBigWig') or 'bedGraphToBigWig'
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("iCLIP_bedtools", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start iCLIP bedtools for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.bed))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"iCLIP_bedtools_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"sample_name=$(basename {output.bed} .bed)\n")
                f.write(f"outdir=$(dirname {output.bed})\n")
                f.write(f"echo \"Processing sample: ${{sample_name}}; output directory: ${{outdir}}\" > {log}\n")
                f.write(f"{params.bedtools} bamtobed -i {input.bam} > {output.bed} 2>> {log}\n")
                f.write(f"\n")
                f.write(f"{params.bedtools} shift -m 1 -p -1 -i {output.bed} -g {input.chromosome_sizes} > ${{outdir}}/${{sample_name}}.shifted.bed 2>> {log}\n")
                f.write(f"{params.bedtools} genomecov -bg -strand + -5 -scale 1000000 -i ${{outdir}}/${{sample_name}}.shifted.bed -g {input.chromosome_sizes} > ${{outdir}}/${{sample_name}}.plus.bedgraph 2>> {log}\n")
                f.write(f"{params.bedtools} genomecov -bg -strand - -5 -scale 1000000 -i ${{outdir}}/${{sample_name}}.shifted.bed -g {input.chromosome_sizes} > ${{outdir}}/${{sample_name}}.minus.bedgraph 2>> {log}\n")
                f.write(f"export LC_COLLATE=C\n")
                f.write(f"sort -k1,1 -k2,2n ${{outdir}}/${{sample_name}}.plus.bedgraph > ${{outdir}}/${{sample_name}}.plus.sorted.bedgraph 2>> {log}\n")
                f.write(f"sort -k1,1 -k2,2n ${{outdir}}/${{sample_name}}.minus.bedgraph > ${{outdir}}/${{sample_name}}.minus.sorted.bedgraph 2>> {log}\n")
                f.write(f"{params.bedGraphToBigWig} ${{outdir}}/${{sample_name}}.plus.sorted.bedgraph {input.chromosome_sizes} ${{outdir}}/${{sample_name}}.plus.bw 2>> {log}\n")
                f.write(f"{params.bedGraphToBigWig} ${{outdir}}/${{sample_name}}.minus.sorted.bedgraph {input.chromosome_sizes} ${{outdir}}/${{sample_name}}.minus.bw 2>> {log}\n")
                f.write(f"rm ${{outdir}}/${{sample_name}}.shifted.bed ${{outdir}}/${{sample_name}}.plus.bedgraph ${{outdir}}/${{sample_name}}.minus.bedgraph ${{outdir}}/${{sample_name}}.plus.sorted.bedgraph ${{outdir}}/${{sample_name}}.minus.sorted.bedgraph\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during iCLIP bedtools for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during iCLIP bedtools for sample {wildcards.sample_id}: {e}")
            raise e
