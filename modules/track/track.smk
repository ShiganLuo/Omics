from snakemake.logging import logger
indir = config.get('indir', "input")
outdir = config.get('outdir', "output")
logdir = config.get('logdir', "log")
samples = config.get('samples', [])
ROOT_DIR = config.get('ROOT_DIR', ".")
rule ucsc_track:
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
        python {params.track_script} -i {input.bigwigs} -o {output.track} > {log} 2>&1
        """
