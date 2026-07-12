include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
peaks_indir = config.get("peaks_indir", "peaks")

samples = config.get("samples", [])

rule frip_score:
    """
    Calculate Fraction of Reads in Peaks (FRiP) for ChIP-seq samples.
    
    FRiP = reads_in_peaks / total_mapped_reads
    
    Uses bedtools intersect to count reads overlapping peaks,
    and samtools flagstat to get total mapped reads.
    """
    input:
        bam = indir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        peak = peaks_indir + "/{sample_id}/{sample_id}_peaks.narrowPeak"
    output:
        frip = outdir + "/{sample_id}/{sample_id}.FRiP.txt"
    log:
        logdir + "/{sample_id}/frip_score.log"
    params:
        bedtools = config.get("Procedure", {}).get("bedtools") or "bedtools",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("frip_score", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start FRiP score calculation for sample {wildcards.sample_id} at {current_time}")
            
            sample_outdir = os.path.dirname(str(output.frip))
            script = os.path.join(sample_outdir, f"frip_score_{current_time}.sh")
            
            script_content = f"""#!/bin/bash
set -euo pipefail

BAM="{input.bam}"
PEAK="{input.peak}"
OUT="{output.frip}"
SAMPLE_ID="{wildcards.sample_id}"
BEDTOOLS="{params.bedtools}"
SAMTOOLS="{params.samtools}"

# Count reads in peaks using bedtools intersect
READS_IN_PEAKS=$($BEDTOOLS intersect -a $BAM -b $PEAK -u | $SAMTOOLS view -c -)

# Get total mapped reads from samtools flagstat
TOTAL_MAPPED=$($SAMTOOLS flagstat $BAM | grep 'mapped (' | grep -v "primary" | head -1 | awk '{{print $1}}')

# Calculate FRiP score
if [ "$TOTAL_MAPPED" -eq 0 ]; then
    FRIP="0.000000"
else
    FRIP=$(awk "BEGIN {{printf \\"%.6f\\", $READS_IN_PEAKS / $TOTAL_MAPPED}}")
fi

# Write output: sample_id\\tfrip_score
echo -e "${{SAMPLE_ID}}\\t${{FRIP}}" > $OUT

echo "FRiP score calculation for sample {wildcards.sample_id} successfully completed!"
echo "Reads in peaks: $READS_IN_PEAKS"
echo "Total mapped reads: $TOTAL_MAPPED"
echo "FRiP score: $FRIP"
"""
            with open(script, "w") as f:
                f.write(script_content)
            
            rule_logger.info(f"Executing script: {script}")
            shell(f"bash {script} > {log_path} 2>&1")
            
            rule_logger.info(f"FRiP score calculation for sample {wildcards.sample_id} completed successfully")
            
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during FRiP score calculation for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during FRiP score calculation for sample {wildcards.sample_id}: {e}")
            raise e


rule frip_score_result:
    """
    Result aggregation rule for subworkflow use rule import.
    """
    input:
        frip = outdir + "/{sample_id}/{sample_id}.FRiP.txt"
