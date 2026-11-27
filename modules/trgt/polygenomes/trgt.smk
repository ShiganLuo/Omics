include: "../../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
karyotype = config.get("Params", {}).get("trgt", {}).get("karyotype") or "XX"
bam_substring = config.get("bam_substring") or ""

_formats_trgt_plot = config.get("Params", {}).get("trgt", {}).get("trgt_plot", {}).get("formats") or ["png"]
if isinstance(_formats_trgt_plot, str):
    _formats_trgt_plot = [_formats_trgt_plot]
_formats_trgt_result = config.get("Params", {}).get("trgt", {}).get("trgt_plot", {}).get("formats") or ["png"]
if isinstance(_formats_trgt_result, str):
    _formats_trgt_result = [_formats_trgt_result]

def get_input_for_trgt_genotype(wildcards):
    logger.info(f"trgt_genotype called with {wildcards}")
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bam.bai")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam.bai")
    in_dict["fasta"] = config["genome"][wildcards.genome]["fasta"]
    in_dict["fai"] = config["genome"][wildcards.genome]["fai"]
    in_dict["bed"] = config["genome"][wildcards.genome]["repeat_bed"]
    return in_dict
rule trgt_genotype:
    input:
        unpack(get_input_for_trgt_genotype)
    output:
        vcf = outdir + "/{genome}/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        bam = outdir + "/{genome}/genotype/{sample_id}/{sample_id}.trgt.spanning.sorted.bam"
    log:
        logdir + "/{genome}/{sample_id}/trgt_genotype.log"
    threads: rule_res("trgt_genotype").get("threads", 4)
    benchmark: logdir + "/{genome}/{sample_id}/trgt_genotype.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "trgt_genotype", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 32000),
    params:
        trgt = config.get("Procedure", {}).get("trgt") or "trgt",
        karyotype = karyotype,
        prefix = lambda wildcards: outdir + f"/{wildcards.genome}/genotype/{wildcards.sample_id}/{wildcards.sample_id}.trgt"
    conda:
        "../trgt.yaml"
    container:
        sif("../trgt.yaml", "trgt_genotype")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start trgt genotype for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
        script = receipt_script(f"trgt_genotype_{wildcards.genome}_{current_time}.sh")
        os.makedirs(os.path.dirname(script), exist_ok=True)
        plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.sample_id}/trgt_genotype.sh")
        cmd = [params.trgt, "genotype",
            "--num-threads", str(threads),
            "--genome", input.fasta,
            "--repeats", input.bed,
            "--karyotype", params.karyotype,
            "--reads", input.bam,
            "--output-prefix", params.prefix
        ]
        with open(script, "w") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(" ".join(cmd) + "\n")
            f.write(f'echo "script completed successfully"\n')
        shutil.copyfile(script, plain_script)
        if SCRIPT_ONLY:
            script_only_touch(output)
        else:
            shell(f"bash {plain_script} >> {log} 2>&1")
rule trgt_plot:
    input:
        vcf = outdir + "/{genome}/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        bam = outdir + "/{genome}/genotype/{sample_id}/{sample_id}.trgt.spanning.sorted.bam",
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"],
        fai = lambda wildcards: config["genome"][wildcards.genome]["fai"],
        bed = lambda wildcards: config["genome"][wildcards.genome]["repeat_bed"]
    output:
        plots = [outdir + "/{genome}/plot/{sample_id}/{sample_id}.trgt.repeat." + fmt for fmt in _formats_trgt_plot],
    log:
        logdir + "/{genome}/{sample_id}/trgt_plot.log"
    benchmark: logdir + "/{genome}/{sample_id}/trgt_plot.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "trgt_plot", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 16000),
    params:
        trgt = config.get("Procedure", {}).get("trgt") or "trgt",
        repeat_id = config.get("Params", {}).get("trgt", {}).get("repeat_id") or "HTT",
        formats = _formats_trgt_plot,
        out_base = outdir + "/{genome}/plot/{sample_id}/{sample_id}.trgt.repeat",
    conda:
        "../trgt.yaml"
    container:
        sif("../trgt.yaml", "trgt_plot")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start trgt plot for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
        script = receipt_script(f"trgt_plot_{wildcards.genome}_{current_time}.sh")
        os.makedirs(os.path.dirname(script), exist_ok=True)
        plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.sample_id}/trgt_plot.sh")
        cmd = [
            params.trgt, "plot",
            "--genome", input.fasta,
            "--repeats", input.bed,
            "--vcf", input.vcf,
            "--spanning-reads", input.bam,
            "--repeat-id", params.repeat_id,
        ]
        with open(script, "w") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            for fmt in params.formats:
                f.write(" ".join(cmd + ["--image", params.out_base + "." + str(fmt)]) + "\n")
            f.write(f'echo "script completed successfully"\n')
        shutil.copyfile(script, plain_script)
        if SCRIPT_ONLY:
            script_only_touch(output)
        else:
            shell(f"bash {plain_script} >> {log} 2>&1")
rule trgt_result:
    input:
        vcf = outdir + "/{genome}/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        plots = [outdir + "/{genome}/plot/{sample_id}/{sample_id}.trgt.repeat." + fmt for fmt in _formats_trgt_result],
