include: "../../../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule bam_sort:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.bam",
        bai = outdir + "/{genome}/{sample_id}/{sample_id}.bam.bai"
    log:
        logdir + "/{genome}/{sample_id}/samtools_sort.log"
    threads: rule_res("bam_sort").get("threads", 8)
    benchmark: logdir + "/{genome}/{sample_id}/samtools_sort.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "bam_sort", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../../samtools.yaml"
    container:
        sif("../../samtools.yaml", "bam_sort")
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools",
        fasta = lambda wildcards: config.get("genome", {}).get(wildcards.genome, {}).get("fasta")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start samtools sort for sample {wildcards.sample_id} at {current_time}")
        cmd1 = [
            params.samtools,
            "sort",
            "-@",
            str(threads),
            "-o",
            output.bam,
            input.bam
        ]
        cmd2 = [
            params.samtools,
            "index",
            "-@",
            str(threads),
            output.bam,
            output.bai
        ]
        script = receipt_script(f"samtools_sort_{wildcards.genome}_{current_time}.sh")
        os.makedirs(os.path.dirname(script), exist_ok=True)
        plain_script = os.path.join(outdir, f"samtools_sort.sh")
        with open(script, "w") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(f'echo "script completed successfully"\n')
        shutil.copyfile(script, plain_script)
        if SCRIPT_ONLY:
            script_only_touch(output)
        else:
            shell("bash {plain_script} >> {log} 2>&1")
rule samtools_sort_index_result:
    input:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.bam",
        bai = outdir + "/{genome}/{sample_id}/{sample_id}.bam.bai"
