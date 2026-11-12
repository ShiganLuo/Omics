include: "../common/common.smk"

# DRclass: DNA(WES somatic) x RNA-seq 联合证据分类模块.
# 所有输入必须为具体文件路径(config: samples.<id>.xxx 显式声明, 无目录扫描).

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

SAMPLES = config.get("samples", {})
PROC = config.get("Procedure", {})
PARAMS = config.get("Params", {}).get("drclass", {})

# genome 兼容两种结构: 扁平 {fasta:...} (模块直跑) 或嵌套
# {default:..., references:{<g>:{fasta:...}}} (subworkflow 传入)
GENOME = config.get("genome", {})
if "fasta" not in GENOME and "references" in GENOME:
    _g = GENOME.get("default") or sorted(GENOME["references"].keys())[0]
    GENOME = GENOME["references"].get(_g, {})


def _s(sid):
    return SAMPLES[sid]


rule drclass_anchors:
    """M4 锚表: somatic PASS(DNA锚) ∪ RNA VCF 检出(RNA锚, 排除DNA锚)."""
    input:
        somatic=lambda w: _s(w.sample_id)["somatic_vcf"],
        rna=lambda w: _s(w.sample_id)["rna_vcf"],
    output:
        tsv=outdir + "/{sample_id}/anchors.tsv",
        bed=outdir + "/{sample_id}/anchors.bed",
    log:
        logdir + "/{sample_id}/drclass_anchors.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=PROC.get("python3") or "python3",
    resources:
        mem_mb=rule_res("drclass_anchors").get("mem_mb", 2000),
    shell:
        r"""
        {params.py_proc} {params.py} anchors \
            --somatic-vcf {input.somatic} \
            --rna-vcf {input.rna} \
            --out-tsv {output.tsv} \
            --out-bed {output.bed} > {log} 2>&1
        """


rule drclass_pileup:
    """M3 RNA 定量: bcftools mpileup 对 RNA BAM 在锚位点计数 (-L anchors.bed)."""
    input:
        bam=lambda w: _s(w.sample_id)["rna_bam"],
        bai=lambda w: _s(w.sample_id)["rna_bam"].replace(".bam", ".bai"),
        fasta=GENOME["fasta"],
        fai=GENOME["fai_index"],
        bed=outdir + "/{sample_id}/anchors.bed",
    output:
        tsv=outdir + "/{sample_id}/rna.pileup.tsv",
    log:
        logdir + "/{sample_id}/drclass_pileup.log"
    params:
        bcftools=PROC.get("bcftools") or "bcftools",
        min_mq=PARAMS.get("min_mq", 20),
        min_bq=PARAMS.get("min_bq", 30),
    resources:
        mem_mb=rule_res("drclass_pileup").get("mem_mb", 8000),
    threads: rule_res("drclass_pileup").get("threads", 4)
    shell:
        r"""
        {params.bcftools} mpileup -f {input.fasta} -R {input.bed} -d 1000 \
            --min-MQ {params.min_mq} --min-BQ {params.min_bq} \
            -a FORMAT/AD,FORMAT/DP -Ou -q 0 -Q 0 {input.bam} 2>> {log} \
        | {params.bcftools} query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n' \
            > {output.tsv} 2>> {log}
        """


rule drclass_classify:
    """状态机判定 -> O1 dna_anchored.tsv + O2 rna_anchored.tsv."""
    input:
        anchors=outdir + "/{sample_id}/anchors.tsv",
        pileup=outdir + "/{sample_id}/rna.pileup.tsv",
        somatic=lambda w: _s(w.sample_id)["somatic_vcf"],
        rna=lambda w: _s(w.sample_id)["rna_vcf"],
        germline=lambda w: _s(w.sample_id).get("germline_vcf"),
        strain=lambda w: _s(w.sample_id).get("strain_vcf", []),
        splice=lambda w: _s(w.sample_id).get("splice_bed", []),
        editing=lambda w: _s(w.sample_id).get("editing_bed", []),
    output:
        dna=outdir + "/{sample_id}/{sample_id}.dna_anchored.tsv",
        rna=outdir + "/{sample_id}/{sample_id}.rna_anchored.tsv",
    log:
        logdir + "/{sample_id}/drclass_classify.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=PROC.get("python3") or "python3",
        min_af=PARAMS.get("min_af", 0.02),
        min_ad=PARAMS.get("min_ad", 4),
        min_dp=PARAMS.get("min_dp", 20),
        germline_af=PARAMS.get("germline_af", 0.20),
        germline_dp=PARAMS.get("germline_dp", 10),
        splice_window=PARAMS.get("splice_window", 5),
    resources:
        mem_mb=rule_res("drclass_classify").get("mem_mb", 4000),
    shell:
        r"""
        {params.py_proc} {params.py} classify \
            --sample {wildcards.sample_id} \
            --anchors {input.anchors} \
            --pileup {input.pileup} \
            --somatic-vcf {input.somatic} \
            --rna-vcf {input.rna} \
            $(if [ -n "{input.germline}" ]; then echo --germline-vcf {input.germline}; fi) \
            $(for f in {input.strain}; do echo -n " --strain-vcf $f"; done) \
            $(for f in {input.splice}; do echo -n " --splice-bed $f"; done) \
            $(for f in {input.editing}; do echo -n " --editing-bed $f"; done) \
            --min-af {params.min_af} --min-ad {params.min_ad} --min-dp {params.min_dp} \
            --germline-af {params.germline_af} --germline-dp {params.germline_dp} \
            --splice-window {params.splice_window} \
            --out-dna {output.dna} --out-rna {output.rna} > {log} 2>&1
        """


rule drclass_summary:
    """O3: 每样本标签计数 + 穿透率 + R4 质控, 全样本汇总."""
    input:
        tables=expand(outdir + "/{s}/{s}.dna_anchored.tsv", s=SAMPLES.keys()) +
               expand(outdir + "/{s}/{s}.rna_anchored.tsv", s=SAMPLES.keys()),
    output:
        tsv=outdir + "/DRclass_summary.tsv",
    log:
        logdir + "/drclass_summary.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=PROC.get("python3") or "python3",
        max_r4=PARAMS.get("max_r4_fraction", 0.30),
    resources:
        mem_mb=rule_res("drclass_summary").get("mem_mb", 2000),
    shell:
        r"""
        {params.py_proc} {params.py} summary \
            --tables {input.tables} \
            --max-r4 {params.max_r4} \
            --out {output.tsv} > {log} 2>&1
        """
