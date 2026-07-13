include: "../common/common.smk"

from snakemake.logging import logger

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

# IGV config shared by igv rules
igv_config = config.get("igv", {})
if not igv_config:
    logger.warning("IGV configuration is missing in config under 'igv'. IGV track rules will fail if triggered.")


rule ucsc_track_single:
    input:
        bigwigs = expand("{indir}/{sample_id}/{sample_id}.bigwig", indir=indir, sample_id=samples)
    output:
        track = outdir + "/ucsc_track.txt"
    log:
        logdir + "/ucsc_track.log"
    conda: "track.yaml"
    params:
        track_script = ROOT_DIR + "/modules/track/bin/track.py"
    run:
        try:
            open(log[0], "w").close()
            rule_logger = setup_logger(logger_name="ucsc_track_bedtools", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ucsc_track_bedtools at {current_time}")

            sample_outdir = os.path.dirname(str(output.track))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"ucsc_track_bedtools_{current_time}.sh")

            cmd = [
                "python", params.track_script,
                "--mode", "ucsc",
                "-i",
            ]
            cmd += list(input.bigwigs)
            cmd += ["-o", str(output.track)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")

            rule_logger.info(f"ucsc_track_bedtools completed successfully at {time.strftime('%Y%m%d_%H%M%S', time.localtime())}")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"ucsc_track_bedtools failed with error: {e}\n")
            logger.error(f"ucsc_track_bedtools failed with error: {e}")
            raise e


rule ucsc_track_iclip:
    input:
        plus_bw = expand("{indir}/{sample_id}/{sample_id}.plus.bw", indir=indir, sample_id=samples),
        minus_bw = expand("{indir}/{sample_id}/{sample_id}.minus.bw", indir=indir, sample_id=samples)
    output:
        track = outdir + "/ucsc_track_iclip.txt"
    log:
        logdir + "/all/ucsc_track_iclip.log"
    conda: "track.yaml"
    params:
        track_script = ROOT_DIR + "/modules/track/bin/track.py"
    run:
        try:
            open(log[0], "w").close()
            rule_logger = setup_logger(logger_name="ucsc_track_iclip", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ucsc_track_iclip at {current_time}")

            sample_outdir = os.path.dirname(str(output.track))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"ucsc_track_iclip_{current_time}.sh")

            cmd = [
                "python", params.track_script,
                "--mode", "ucsc",
                "-i",
            ]
            cmd += list(input.plus_bw)
            cmd += list(input.minus_bw)
            cmd += ["-o", str(output.track)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")

            rule_logger.info(f"ucsc_track_iclip completed successfully at {time.strftime('%Y%m%d_%H%M%S', time.localtime())}")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"ucsc_track_iclip failed with error: {e}\n")
            logger.error(f"ucsc_track_iclip failed with error: {e}")
            raise e


rule igv_track_single:
    input:
        bigwigs = expand("{indir}/{sample_id}/{sample_id}.bigwig", indir=indir, sample_id=samples)
    output:
        html = outdir + "/igv_track.html"
    log:
        logdir + "/all/igv_track.log"
    conda: "track.yaml"
    params:
        track_script = ROOT_DIR + "/modules/track/bin/track.py",
        igv_config = igv_config,
        config_json = outdir + "/igv_track_bedtools_config.json",
    run:
        import json
        if not params.igv_config:
            raise ValueError(
                "IGV configuration is missing in config under 'igv'. "
                "Please provide the necessary configuration for IGV track generation."
            )
        try:
            open(log[0], "w").close()
            rule_logger = setup_logger(logger_name="igv_track_bedtools", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start igv_track_bedtools at {current_time}")

            sample_outdir = os.path.dirname(str(output.html))
            os.makedirs(sample_outdir, exist_ok=True)
            config_json_path = str(params.config_json)
            script = os.path.join(sample_outdir, f"igv_track_bedtools_{current_time}.sh")

            rule_logger.info(f"Writing IGV config JSON to {config_json_path}")
            with open(config_json_path, "w") as f:
                json.dump(params.igv_config, f, indent=2)

            cmd = [
                "python", params.track_script,
                "--mode", "igv",
                "--config", config_json_path,
                "-i",
            ]
            cmd += list(input.bigwigs)
            cmd += ["-o", str(output.html)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")

            rule_logger.info(f"igv_track_bedtools completed successfully at {time.strftime('%Y%m%d_%H%M%S', time.localtime())}")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"igv_track_bedtools failed with error: {e}\n")
            logger.error(f"igv_track_bedtools failed with error: {e}")
            raise e


rule igv_track_iclip:
    input:
        plus_bw = expand("{indir}/{sample_id}/{sample_id}.plus.bw", indir=indir, sample_id=samples),
        minus_bw = expand("{indir}/{sample_id}/{sample_id}.minus.bw", indir=indir, sample_id=samples)
    output:
        html = outdir + "/igv_track_iclip.html"
    log:
        logdir + "/all/igv_track_iclip.log"
    conda: "track.yaml"
    params:
        track_script = ROOT_DIR + "/modules/track/bin/track.py",
        igv_config = igv_config,
        config_json = outdir + "/igv_track_iclip_config.json",
    run:
        import json
        if not params.igv_config:
            raise ValueError(
                "IGV configuration is missing in config under 'igv'. "
                "Please provide the necessary configuration for IGV track generation."
            )
        try:
            open(log[0], "w").close()
            rule_logger = setup_logger(logger_name="igv_track_iclip", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start igv_track_iclip at {current_time}")

            sample_outdir = os.path.dirname(str(output.html))
            os.makedirs(sample_outdir, exist_ok=True)
            config_json_path = str(params.config_json)
            script = os.path.join(sample_outdir, f"igv_track_iclip_{current_time}.sh")

            rule_logger.info(f"Writing IGV config JSON to {config_json_path}")
            with open(config_json_path, "w") as f:
                json.dump(params.igv_config, f, indent=2)

            cmd = [
                "python", params.track_script,
                "--mode", "igv",
                "--config", config_json_path,
                "-i",
            ]
            cmd += list(input.plus_bw)
            cmd += list(input.minus_bw)
            cmd += ["-o", str(output.html)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")

            rule_logger.info(f"igv_track_iclip completed successfully at {time.strftime('%Y%m%d_%H%M%S', time.localtime())}")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"igv_track_iclip failed with error: {e}\n")
            logger.error(f"igv_track_iclip failed with error: {e}")
            raise e
