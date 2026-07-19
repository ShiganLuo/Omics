include: "../common/common.smk"

MODULE_DIR = os.path.join(ROOT_DIR, "modules", "annovar")
logger.info(f"Annovar module directory: {MODULE_DIR}")

# need test
rule TEcoutCPM:
    input:
       infile = outdir + "/counts/humanTEcount.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCPM.cntTable"
    log:
        log = outdir + "/log/human/TEcoutCPM.log"
    params:
        script = MODULE_DIR + "/scripts/SNP/run-NormCountMat.R",
        Rscript = config["Procedure"]["Rscript"]
    conda:
        "annovar.yaml"
    run:
        log_path = str(log.log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="TEcoutCPM_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start TEcoutCPM run at {current_time}")
            script = os.path.join(outdir, f"TEcoutCPM_{current_time}.sh")
            cmd = [params.Rscript, params.script, "-i", input.infile, "-o", output.outfile]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")

rule commonExpression:
    input:
        infile = outdir + "/counts/humanTEcountCPM.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCommon.cntTable"
    log:
        log = outdir + "/log/human/commonExpression.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = MODULE_DIR + "scripts/SNP/commonExpression.py"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="commonExpression", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start commonExpression at {current_time}")
            script = os.path.join(outdir, f"commonExpression_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--input", input.infile,
                "--output", output.outfile,
                "--threshold", "5"
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"commonExpression failed with error: {e}\n")
            raise f"commonExpression failed with error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"commonExpression completed at {current_time}")

rule getBed:
    input:
        infile = outdir + "/counts/humanTEcountCommon.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCommon.bed"
    log:
        log = outdir + "/log/human/getBed.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = MODULE_DIR + "scripts/SNP/getBed.py",
        gtf = config['bed']['human']['gtf'],
        TE_gtf = config['bed']['human']['TE_gtf']
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="getBed", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start getBed at {current_time}")
            script = os.path.join(outdir, f"getBed_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--input", input.infile,
                "--output", output.outfile,
                "--Gtf", params.gtf,
                "--TEGtf", params.TE_gtf
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"getBed failed with error: {e}\n")
            raise f"getBed failed with error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"getBed completed at {current_time}")

rule vcfIntersectBed:
    input:
        vcf = outdir + "/filter/vcf/human/{sample_id}.vcf.gz",
        bed = outdir + "/counts/humanTEcountCommon.bed"
    output:
        outfile = outdir + "/filter/vcf/human/{sample_id}Common.vcf"
    log:
        log = outdir + "/log/human/{sample_id}/vcfIntersectBed.log"
    conda:
        config['conda']['RNA-SNP']
    threads:4 #防止同时执行太多，爆内存
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="vcfIntersectBed", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start vcfIntersectBed for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/vcfIntersectBed_{current_time}.sh")
            cmd = [
                "bedtools", "intersect",
                "-a", input.vcf,
                "-b", input.bed,
                "-wa", "-wb"
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"vcfIntersectBed failed for sample {wildcards.sample_id} with error: {e}\n")
            raise f"vcfIntersectBed failed for sample {wildcards.sample_id} with error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"vcfIntersectBed completed for sample {wildcards.sample_id} at {current_time}")

rule annovar_convert:
    input:
        vcf = outdir + "/filter/vcf/human/{sample_id}Common.vcf"
    output:
        avinput = outdir + "/annovar/human/{sample_id}/{sample_id}.avinput"
    log:
        log = outdir + "/log/human/{sample_id}/annovar_convert.log"
    params:
        convert = "/opt/annovar/convert2annovar.pl",
        perl = config["Procedure"]["perl"]
    threads:4 #防止同时执行太多，爆内存
    conda:
        "annovar.yaml"
    run:
        log_path = str(log.log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="annovar_convert_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start annovar_convert run for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/annovar_convert_{current_time}.sh")
            cmd = [params.perl, params.convert, "-format", "vcf4", "-withfreq", input.vcf, ">", output.avinput, "2>", log_path]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")

rule annovar_table:
    input:
        avinput = outdir + "/annovar/human/{sample_id}/{sample_id}.avinput"
    output:
        outfile = outdir + "/annovar/human/{sample_id}/{sample_id}.GRCh38_multianno.csv"
    log:
        log = outdir + "/log/human/{sample_id}/annovar_table.log"
    params:
        db = config['annovar']['human']['db'],
        buildver = config['annovar']['human']['buildver'],
        # annotate = "/opt/annovar/annotate_variation.pl",
        table = config["Procedure"]["table_annovar"],
        out = outdir + "/annovar/human/{sample_id}/{sample_id}"
    threads: 4 #防止同时执行太多，爆内存
    conda:
        "annovar.yaml"
    run:
        log_path = str(log.log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="annovar_table_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start annovar_table run for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/annovar_table_{current_time}.sh")
            cmd = [params.perl, params.table, input.avinput, params.db, "-buildver", params.buildver, "-out", params.out, "-remove", "-protocol", "refGene", "-operation", "g", "-nastring", ".", "-csvout"]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")
