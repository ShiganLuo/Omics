include: "../common/common.smk"
indir = config.get('indir', "input")
outdir = config.get('outdir', "output")
logdir = config.get('logdir', "log")

rule samtools_dedup:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bam = outdir + "/{sample_id}/{sample_id}.dedup.bam",
        bai = outdir + "/{sample_id}/{sample_id}.dedup.bam.bai",
    log:
        logdir + "/{sample_id}/samtools_dedup.log"
    threads: 12
    conda:
        "igv.yaml"
    params:
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    run:
        log_path = str(log)
        try:
            rule_logger = setup_logger("dedup",log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start rule dedup for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(sample_outdir, f"samtools_dedup_{current_time}.sh")
            cmd = [
                params.samtools, "sort",
                "-n",
                "-@", str(threads),
                input.bam, "|",
                params.samtools, "fixmate",
                "-m", "-", "-", "|",
                params.samtools, "sort",
                "-@", str(threads),
                "-", "|",
                params.samtools, "markdup",
                "-r", 
                "-@", str(threads),
                "-",
                output.bam, "&&",
                params.samtools, "index",
                "-@", str(threads),
                output.bam
            ]
            success_echo = f'echo "dedup for sample {wildcards.sample_id} successfully completed !"'
            with open(script_path,"w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(success_echo + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"samtools dedup failed: {e}\n")
            logger.error(f"samtools dedup failed: {e}\n")
            raise f"samtools dedup failed: {e}\n"
            

rule wig:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.dedup.bam",
        bai = outdir + "/{sample_id}/{sample_id}.dedup.bam.bai"
    output:
        bigwig = outdir + "/{sample_id}/{sample_id}.bigwig"
    log:
        log = logdir + "/{sample_id}/wig.log"
    conda:
        "igv.yaml"
    threads: 12 
    params:
        binSize= config.get('Params',{}).get('bamCoverage',{}).get('binSize') or 50,
        bamCoverage = config.get('Procedure',{}).get('bamCoverage') or 'bamCoverage',
        normalizeUsing = config.get('Params', {}).get('bamCoverage',{}).get('normalizeUsing') or "CPM",
        offset = config.get('Params', {}).get('bamCoverage',{}).get('offset') or None,
        extendReads = config.get('Params', {}).get('bamCoverage',{}).get('extendReads') or False
    run:
        log_path = str(log)
        try:
            rule_logger = setup_logger("wig",log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start rule wig for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(sample_outdir, f"bamCoverage_wig_{current_time}.sh")
            cmd = [
                params.bamCoverage,
                "--numberOfProcessors", str(threads),
                "--binSize", str(params.binSize),
                "--normalizeUsing", params.normalizeUsing,
                "-b", input.bam,
                "-o", output.bigwig
            ]
            if isinstance(params.extendReads, bool):
                if params.extendReads:
                    cmd += ["--extendReads"]
            elif isinstance(params.extendReads, int) and params.extendReads > 0:
                cmd += ["--extendReads", str(params.extendReads)]
            if params.offset:
                cmd += ["--Offset", params.offset]
            with open(script_path,"w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"bamCoverage wig failed: {e}\n")
            logger.error(f"bamCoverage wig failed: {e}\n")
            raise f"bamCoverage wig failed: {e}\n"

rule igv_result:
    input:
        bigwig = outdir + "/{sample_id}/{sample_id}.bigwig"
