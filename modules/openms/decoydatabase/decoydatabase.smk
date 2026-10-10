include: "../../common/common.smk"

indir = config.get("indir", "data/fasta")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
protein_fasta = config.get("genome", {}).get("protein_fasta")
if not protein_fasta or not os.path.exists(protein_fasta):
    raise ValueError("Protein FASTA file path is not specified in the configuration under 'genome.protein_fasta' or does not exist.")

# Get parameters
decoy_string = config.get("Params", {}).get("decoy_database", {}).get("decoy_string", "DECOY_")
decoy_string_position = config.get("Params", {}).get("decoy_database", {}).get("decoy_string_position", "prefix")
method = config.get("Params", {}).get("decoy_database", {}).get("method", "shuffle")
shuffle_max_attempts = config.get("Params", {}).get("decoy_database", {}).get("shuffle_max_attempts", 30)
shuffle_sequence_identity_threshold = config.get("Params", {}).get("decoy_database", {}).get("shuffle_sequence_identity_threshold", 0.5)

# Get OpenMS executable
openms = config.get("Procedure", {}).get("openms") or "DecoyDatabase"

rule decoy_database:
    input:
        protein_fasta = protein_fasta
    output:
        decoy_fasta = outdir + "/genome_decoy.fasta"
    log:
        logdir + "/decoy_database/decoy_database.log"
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        openms = openms,
        decoy_string = decoy_string,
        decoy_string_position = decoy_string_position,
        method = method,
        shuffle_max_attempts = shuffle_max_attempts,
        shuffle_sequence_identity_threshold = shuffle_sequence_identity_threshold
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("decoy_database",log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start decoy database generation at {current_time}")
            # Guard: a nucleic-acid FASTA here silently yields zero PSMs later.
            seen = set()
            with open(input.protein_fasta) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if seen:
                            break
                        continue
                    seen.update(line.strip().upper())
                    if len(seen) > 5:
                        break
            if seen and seen <= set("ACGTN"):
                raise ValueError(
                    f"{input.protein_fasta} looks like a NUCLEIC-ACID FASTA (sequence "
                    f"alphabet is only ACGTN). Protein search needs a protein "
                    f"FASTA (e.g. GENCODE pep.all / UniProt proteome). Set "
                    f"'genome.protein_fasta' to the protein database."
                )
            script = os.path.join(outdir, f"decoy_database_{current_time}.sh")
            cmd = [
                params.openms,
                "-in", input.protein_fasta,
                "-out", output.decoy_fasta,
                "-decoy_string", params.decoy_string,
                "-decoy_string_position", params.decoy_string_position,
                "-method", params.method,
                "-shuffle_max_attempts", str(params.shuffle_max_attempts),
                "-shuffle_sequence_identity_threshold", str(params.shuffle_sequence_identity_threshold)
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'Decoy database generation for {input.protein_fasta} was completed successfully'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule decoy_database was call failed,error: {e}")
            raise RuntimeError(f"rule decoy_database was call failed,error: {e}")

rule decoy_database_result:
    input:
        decoy_fasta = outdir + "/genome_decoy.fasta"
