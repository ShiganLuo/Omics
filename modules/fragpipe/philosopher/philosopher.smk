include: "../../common/common.smk"

indir = config.get("indir", "data/msfragger")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
search_fasta = config.get("genome", {}).get("search_fasta", "")

# Get parameters
phil_params = config.get("Params", {}).get("philosopher", {})
pep_fdr = phil_params.get("pep_fdr", 0.01)
prot_fdr = phil_params.get("prot_fdr", 0.01)
ion_fdr = phil_params.get("ion_fdr", 0.01)
picked = phil_params.get("picked", True)

# Get executables
philosopher = config.get("Procedure", {}).get("philosopher") or "philosopher"

rule philosopher_database:
    input:
        fasta = lambda wildcards: _require(search_fasta, "genome.search_fasta", "philosopher"),
    output:
        db = outdir + "/database/fragpipe_db.fasta",
    log:
        logdir + "/philosopher_database/philosopher_database.log"
    threads: 1
    conda:
        "../fragpipe.yaml"
    container:
        sif("../fragpipe.yaml")
    params:
        philosopher = philosopher,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("philosopher_database", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start Philosopher database build at {current_time}")
            db_dir = os.path.join(outdir, "database")
            os.makedirs(db_dir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"philosopher_database_{current_time}.sh")
            plain_script = os.path.join(db_dir, f"philosopher_database.sh")
            # Philosopher writes a dated '<date>-decoys-contam-<name>.fas';
            # rename to a stable output so reruns stay mtime-stable.
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"cd {shlex.quote(db_dir)}\n")
                f.write(f"{shlex.quote(str(params.philosopher))} workspace --init\n")
                f.write(f"{shlex.quote(str(params.philosopher))} database --custom {shlex.quote(str(input.fasta))} "
                        "--contam --prefix rev_ --decoymode 0\n")
                f.write("gen=$(ls -1 *-decoys-contam-*.fas | head -1)\n")
                f.write(f"mv \"$gen\" {shlex.quote(str(output.db))}\n")
                f.write(f"echo 'Philosopher database build completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during Philosopher database build: {str(e)}\n")
            raise RuntimeError(f"Error during Philosopher database build: {str(e)}\n")

rule philosopher_filter:
    input:
        pepxml = indir + "/{sample_id}/{sample_id}.pepXML",
        db = outdir + "/database/fragpipe_db.fasta",
    output:
        validated = outdir + "/{sample_id}/interact-{sample_id}.pep.xml",
        psm = outdir + "/{sample_id}/psm.tsv",
        peptide = outdir + "/{sample_id}/peptide.tsv",
    log:
        logdir + "/{sample_id}/philosopher_filter.log"
    threads: 1
    conda:
        "../fragpipe.yaml"
    container:
        sif("../fragpipe.yaml")
    params:
        philosopher = philosopher,
        pep_fdr = pep_fdr,
        prot_fdr = prot_fdr,
        ion_fdr = ion_fdr,
        picked = picked,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("philosopher_filter", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start Philosopher validation/filter for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            local_pepxml = os.path.join(sample_outdir, f"{wildcards.sample_id}.pepXML")
            if not os.path.exists(local_pepxml):
                try:
                    os.link(str(input.pepxml), local_pepxml)
                except OSError:
                    shutil.copy(str(input.pepxml), local_pepxml)
            script = os.path.join(RECEIPT_DIR, f"philosopher_filter_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"philosopher_filter.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"cd {shlex.quote(sample_outdir)}\n")
                f.write(f"{shlex.quote(str(params.philosopher))} workspace --init\n")
                # --annotate rebuilds the per-workspace .meta/db.bin the
                # protein-level filtering reads.
                f.write(f"{shlex.quote(str(params.philosopher))} database --annotate {shlex.quote(str(input.db))}\n")
                f.write(f"{shlex.quote(str(params.philosopher))} peptideprophet --database {shlex.quote(str(input.db))} "
                        f"--decoy rev_ --nonparam --expectscore --accmass --decoyprobs "
                        f"{shlex.quote(local_pepxml)}\n")
                filter_cmd = [
                    str(params.philosopher), "filter",
                    "--pepxml", str(output.validated),
                    "--pep", str(params.pep_fdr),
                    "--prot", str(params.prot_fdr),
                    "--ion", str(params.ion_fdr),
                ]
                if params.picked:
                    filter_cmd.append("--picked")
                f.write(" ".join(shlex.quote(str(x)) for x in filter_cmd) + "\n")
                f.write(f"{shlex.quote(str(params.philosopher))} report\n")
                f.write(f"echo 'Philosopher validation/filter for {wildcards.sample_id} completed at {current_time}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during Philosopher filter: {str(e)}\n")
            raise RuntimeError(f"Error during Philosopher filter: {str(e)}\n")

rule philosopher_result:
    input:
        expand(outdir + "/{sample_id}/psm.tsv", sample_id=samples),
        outdir + "/database/fragpipe_db.fasta",
