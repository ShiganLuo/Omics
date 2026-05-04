indir = config.get('indir', "input")
outdir = config.get('outdir', "output")
logdir = config.get('logdir', "log")
samples = config.get('samples', [])
ROOT_DIR = config.get('ROOT_DIR', ".")


rule ucsc_track_bedtools:
    input:
        bigwigs = expand("{indir}/{sample_id}.bigwig", indir=indir, sample_id=samples)
    output:
        track = outdir + "/ucsc_track.txt"
    log:
        logdir + "/ucsc_track.log"
    conda:
        "track.yaml"
    params:
        track_script = ROOT_DIR + "/modules/track/bin/track.py"
    shell:
        """
        python {params.track_script} --mode ucsc -i {input.bigwigs} -o {output.track} > {log} 2>&1
        """

rule ucsc_track_iclip:
    input:
        plus_bw = expand("{indir}/{sample_id}/{sample_id}.plus.bw", indir=indir, sample_id=samples),
        minus_bw = expand("{indir}/{sample_id}/{sample_id}.minus.bw", indir=indir, sample_id=samples)
    output:
        track = outdir + "/ucsc_track_iclip.txt"
    log:
        logdir + "/ucsc_track_iclip.log"
    conda:
        "track.yaml"
    params:
        track_script = ROOT_DIR + "/modules/track/bin/track.py"
    shell:
        """
        python {params.track_script} --mode ucsc -i {input.plus_bw} {input.minus_bw} -o {output.track} > {log} 2>&1
        """

rule igv_track_bedtools:
    input:
        bigwigs = expand("{indir}/{sample_id}.bigwig", indir=indir, sample_id=samples)
    output:
        html = outdir + "/igv_track.html"
    log:
        logdir + "/igv_track.log"
    conda:
        "track.yaml"
    params:
        track_script = ROOT_DIR + "/modules/track/bin/track.py",
        igv_config = config.get('igv', {}),
        config_json = outdir + "/igv_track_config.json",
    run:
        import json
        if not params.igv_config:
            logger.error("IGV configuration is missing in the config file under 'Params: igv'. Please provide the necessary configuration for IGV track generation.")
            raise ValueError("Missing IGV configuration")
        with open(log[0], "w") as lf:
            lf.write("=== IGV track generation start ===\n")

            lf.write("Writing config JSON...\n")
            with open(params.config_json, "w") as f:
                json.dump(params.igv_config, f, indent=2)

            lf.write(f"Config written to: {params.config_json}\n")

        shell(
        """
        python {params.track_script} \
            --mode igv \
            --config {params.config_json} \
            -i {input.plus_bw} {input.minus_bw} \
            -o {output.html} >> {log} 2>&1
        """
        )

rule igv_track_iclip:
    input:
        plus_bw = expand("{indir}/{sample_id}/{sample_id}.plus.bw", indir=indir, sample_id=samples),
        minus_bw = expand("{indir}/{sample_id}/{sample_id}.minus.bw", indir=indir, sample_id=samples)
    output:
        html = outdir + "/igv_track_iclip.html"
    log:
        logdir + "/igv_track_iclip.log"
    conda:
        "track.yaml"
    params:
        track_script = ROOT_DIR + "/modules/track/bin/track.py",
        igv_config = config.get('igv', {}),
        config_json = outdir + "/igv_track_config.json",
    run:
        import json
        if not params.igv_config:
            logger.error("IGV configuration is missing in the config file under 'Params: igv'. Please provide the necessary configuration for IGV track generation.")
            raise ValueError("Missing IGV configuration")
        with open(log[0], "w") as lf:
            lf.write("=== IGV track generation start ===\n")

            lf.write("Writing config JSON...\n")
            with open(params.config_json, "w") as f:
                json.dump(params.igv_config, f, indent=2)

            lf.write(f"Config written to: {params.config_json}\n")

        shell(
        """
        python {params.track_script} \
            --mode igv \
            --config {params.config_json} \
            -i {input.plus_bw} {input.minus_bw} \
            -o {output.html} >> {log} 2>&1
        """
        )
