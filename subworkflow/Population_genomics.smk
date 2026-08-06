"""Orchestrate atomic population-genomics modules."""

from itertools import combinations

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
sample_bams = config.get("sample_bams", {})
populations = config.get("populations", {})
analysis = config.get("Params", {}).get("analysis", {})
pop_names = sorted(populations)

filtered_vcf = f"{outdir}/variants/filtered/population_snps.vcf.gz"
filtered_tbi = f"{filtered_vcf}.tbi"
plink_prefix = f"{outdir}/analysis/plink/population_snps"
pop_pairs = list(combinations(pop_names, 2))
outfiles = config.get("outfiles", [])
if not outfiles:
    outfiles = [filtered_vcf, filtered_tbi]
    if analysis.get("pca", True):
        outfiles += [f"{plink_prefix}.eigenvec", f"{plink_prefix}.eigenval"]
    if analysis.get("structure", True):
        outfiles += [f"{outdir}/analysis/structure/admixture.K{int(k)}.Q" for k in analysis.get("admixture_k", [2, 3, 4, 5, 6])]
    if analysis.get("diversity", True):
        for population in pop_names:
            outfiles += [f"{outdir}/analysis/diversity/{population}.windowed.pi", f"{outdir}/analysis/diversity/{population}.Tajima.D"]
    if analysis.get("fst", True):
        outfiles += [f"{outdir}/analysis/fst/{a}_vs_{b}.weir.fst" for a, b in pop_pairs]
    if analysis.get("history", False):
        outfiles.append(f"{outdir}/analysis/history/easySFS.done")
    if analysis.get("ld_decay", False):
        outfiles += [f"{outdir}/analysis/ld_decay/{population}.stat.gz" for population in pop_names]
    if analysis.get("gwas", False):
        outfiles.append(f"{outdir}/analysis/gwas/gwas.done")

base_config = dict(config)
base_config.update({"ROOT_DIR": ROOT_DIR, "indir": indir, "outdir": outdir, "logdir": logdir, "samples": samples, "sample_bams": sample_bams, "populations": populations, "outfiles": outfiles})

metadata_config = {"ROOT_DIR": ROOT_DIR, "outdir": outdir, "logdir": logdir, "populations": populations}
module population_metadata:
    snakefile: "../modules/population_metadata/population_metadata.smk"
    config: metadata_config
use rule population_metadata from population_metadata as Population_metadata

gatk_population_config = dict(base_config)
module gatk_population:
    snakefile: "../modules/gatk/gatk_population/gatk_population.smk"
    config: gatk_population_config
use rule gatk_population_haplotype_caller from gatk_population as Population_gatk_haplotype_caller
use rule gatk_population_joint_genotyping from gatk_population as Population_gatk_joint_genotyping

bcftools_population_config = dict(base_config)
bcftools_population_config["input_vcf"] = f"{outdir}/variants/joint/joint.vcf.gz"
module bcftools_population:
    snakefile: "../modules/bcftools/bcftools_population/bcftools_population.smk"
    config: bcftools_population_config
use rule bcftools_population_filter from bcftools_population as Population_bcftools_filter

plink2_population_config = dict(base_config)
plink2_population_config["input_vcf"] = filtered_vcf
module plink2_population:
    snakefile: "../modules/plink2/plink2_population.smk"
    config: plink2_population_config
use rule * from plink2_population as Population_plink2_*

vcftools_population_config = dict(base_config)
vcftools_population_config["input_vcf"] = filtered_vcf
module vcftools_population:
    snakefile: "../modules/vcftools/vcftools_population.smk"
    config: vcftools_population_config
use rule * from vcftools_population as Population_vcftools_*

if analysis.get("structure", True):
    admixture_config = dict(base_config)
    module admixture_population:
        snakefile: "../modules/admixture/admixture_population.smk"
        config: admixture_config
    use rule admixture_population from admixture_population as Population_admixture

if analysis.get("history", False):
    easy_sfs_config = dict(base_config)
    easy_sfs_config["input_vcf"] = filtered_vcf
    module easySFS_population:
        snakefile: "../modules/easySFS/easySFS_population.smk"
        config: easy_sfs_config
    use rule easySFS_population_history from easySFS_population as Population_easySFS_history

if analysis.get("ld_decay", False):
    poplddecay_config = dict(base_config)
    poplddecay_config["input_vcf"] = filtered_vcf
    module PopLDdecay_population:
        snakefile: "../modules/PopLDdecay/PopLDdecay_population.smk"
        config: poplddecay_config
    use rule PopLDdecay_population from PopLDdecay_population as Population_PopLDdecay

rule all:
    input:
        outfiles
