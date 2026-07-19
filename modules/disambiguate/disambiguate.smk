include: "../common/common.smk"
from typing import List
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
paired_samples =  config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
genome_pairs: List[str] = config.get("genome_pairs", [])
genomeA, genomeB = genome_pairs


def get_inputFile_for_ngs_disambiguate(wildcards):
    logger.info(f"[get_inputFile_for_ngs_disambiguate] called with wildcards: {wildcards}")
    return {
        "bamA": f"{indir}/{genomeA}/{wildcards.sample_id}/{wildcards.sample_id}.bam",
        "bamB": f"{indir}/{genomeB}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    }

rule ngs_disambiguate:
    input:
        bamA = lambda wc: get_inputFile_for_ngs_disambiguate(wc)["bamA"],
        bamB = lambda wc: get_inputFile_for_ngs_disambiguate(wc)["bamB"]
    output:
        raw_bamA = temp(outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesA.bam"),
        raw_bamB = temp(outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesB.bam"),
        raw_ambiguousA = temp(outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesA.bam"),
        raw_ambiguousB = temp(outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesB.bam"),
        summary = outdir + "/{sample_id}/{sample_id}_summary.tsv"
    params:
        bamA_sortN = temp(outdir + "/{sample_id}/{sample_id}.bamA.sortN.bam"),
        bamB_sortN = temp(outdir + "/{sample_id}/{sample_id}.bamB.sortN.bam"),
        outdir = lambda wc: f"{outdir}/{wc.sample_id}",
        aligner = config.get("Params", {}).get("ngs_disambiguate", {}).get("aligner") or "hisat2",
        ngs_disambiguate = config.get("Procedure", {}).get("ngs_disambiguate") or "ngs_disambiguate",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    threads: 4
    conda:
        "disambiguate.yaml"
    log:
        logdir + "/{sample_id}/ngs_disambiguate.log"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("ngs_disambiguate", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ngs_disambiguate for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(params.outdir, f"ngs_disambiguate_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.samtools} sort -n -@ {threads} -o {params.bamA_sortN} {input.bamA} 2>> {log_path}\n")
                f.write(f"{params.samtools} sort -n -@ {threads} -o {params.bamB_sortN} {input.bamB} 2>> {log_path}\n")
                f.write(f"\n")
                f.write(f"{params.ngs_disambiguate} \\\n")
                f.write(f"    -s {wildcards.sample_id} \\\n")
                f.write(f"    -o {params.outdir} \\\n")
                f.write(f"    -a {params.aligner} \\\n")
                f.write(f"    {params.bamA_sortN} \\\n")
                f.write(f"    {params.bamB_sortN}\n")
                f.write(f"rm -f {params.bamA_sortN} {params.bamB_sortN}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during ngs_disambiguate for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during ngs_disambiguate for sample {wildcards.sample_id}: {e}")
            raise e

rule disambiguate_sort_rename:
    input:
        summary = outdir + "/{sample_id}/{sample_id}_summary.tsv",
        raw_bamA = outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesA.bam",
        raw_bamB = outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesB.bam",
        raw_ambiguousA = outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesA.bam",
        raw_ambiguousB = outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesB.bam"
    output:
        clean_bamA = outdir + "/{sample_id}/{sample_id}" + f".disambiguatedSpecies_{genome_pairs[0]}.bam",
        clean_bamB = outdir + "/{sample_id}/{sample_id}" + f".disambiguatedSpecies_{genome_pairs[1]}.bam",
        ambiguous_bamA = outdir + "/{sample_id}/{sample_id}" + f".ambiguousSpecies_{genome_pairs[0]}.bam",
        ambiguous_bamB = outdir + "/{sample_id}/{sample_id}" + f".ambiguousSpecies_{genome_pairs[1]}.bam",
        clean_summary = outdir + "/{sample_id}/{sample_id}_summary_renamed.tsv"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools",
        speciesA = genome_pairs[0],
        speciesB = genome_pairs[1]
    threads: 4
    conda:
        "disambiguate.yaml"
    log:
        logdir + "/{sample_id}/sort_rename.log"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("disambiguate_sort_rename", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start disambiguate_sort_rename for sample {wildcards.sample_id} at {current_time}")
            outdir_sample = os.path.dirname(str(output.clean_bamA))
            script = os.path.join(outdir_sample, f"sort_rename_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"sed '1s/unique species A pairs/unique species {params.speciesA} pairs/; \\\n")
                f.write(f"    1s/unique species B pairs/unique species {params.speciesB} pairs/' {input.summary} > {output.clean_summary}\n")
                f.write(f"{params.samtools} sort -@ {threads} -o {output.clean_bamA} {input.raw_bamA}\n")
                f.write(f"{params.samtools} sort -@ {threads} -o {output.clean_bamB} {input.raw_bamB}\n")
                f.write(f"\n")
                f.write(f"{params.samtools} sort -@ {threads} -o {output.ambiguous_bamA} {input.raw_ambiguousA}\n")
                f.write(f"{params.samtools} sort -@ {threads} -o {output.ambiguous_bamB} {input.raw_ambiguousB}\n")
                f.write(f"\n")
                f.write(f"{params.samtools} index {output.clean_bamA}\n")
                f.write(f"{params.samtools} index {output.clean_bamB}\n")
                f.write(f"rm -f {input.summary}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during disambiguate_sort_rename for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during disambiguate_sort_rename for sample {wildcards.sample_id}: {e}")
            raise e

rule disambiguate_report:
    input:
        reports = expand(outdir + "/{sample_id}/{sample_id}_summary_renamed.tsv", sample_id=paired_samples + single_samples)
    output:
        report = outdir + "/disambiguate_qc.tsv"
    params:
        combine_script = ROOT_DIR + "/modules/disambiguate/combineDisambiguateQC.py"
    log:
        logdir + "/disambiguate_report.log"
    conda:
        "disambiguate.yaml"
    threads: 1
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("disambiguate_report", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start disambiguate_report at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            script = os.path.join(report_dir, f"disambiguate_report_{current_time}.sh")
            reports_str = " ".join(input.reports)
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python {params.combine_script} \\\n")
                f.write(f"    -i {reports_str} \\\n")
                f.write(f"    -o {output.report}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during disambiguate_report: {e}\n")
            logger.error(f"Error occurred during disambiguate_report: {e}")
            raise e

rule ngs_disambiguate_result:
    input:
        bamA = lambda wc: f"{outdir}/{wc.sample_id}/{wc.sample_id}.disambiguatedSpecies_{genome_pairs[0]}.bam",
        bamB = lambda wc: f"{outdir}/{wc.sample_id}/{wc.sample_id}.disambiguatedSpecies_{genome_pairs[1]}.bam"
