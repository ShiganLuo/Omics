import pandas as pd
import logging
import sys
import re
from pathlib import Path
from collections import defaultdict
from typing import DefaultDict, List, Dict
import subprocess
import os
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,  # output: stdout ,not stderr
	datefmt='%Y-%m-%d %H:%M:%S'
)


class SNPMetadata:
    def __init__(self,meta:str,fqDir:str):
        self.meta = meta
        self.fqDir = Path(fqDir)
    
    def loadMeta(self):
        df = pd.read_csv(self.meta,sep="\t")
        requeired_column = ["data_id","sample_id","organism"]

        missing_cols = [col for col in requeired_column if col not in df.columns]
        if missing_cols:
            logging.error(f"The columns: {missing_cols} must be contained in metadata file")
            logging.info(f"the read columns actually is : {df.columns}")
            raise ValueError(f"The columns: {missing_cols} must be contained in metadata file")
        else:
            logging.info(f"all columns is contained in metadata file: {requeired_column}")
        df = df[requeired_column]
        return df
    
    def isPairedEndSample(self,sample:str) -> str:
        """
        Determine the sequencing strategy type for a specified sample
        PAIRED: 
            _R1*fastq , _R2*fastq or _R1*fq , _R2*fq
            _1*fastq , _2*fastq or _1*fq , _2*fq
        SINGLE: not PAIRED
        """
        fqDir = self.fqDir
        fastq_files = list(fqDir.glob(f"{sample}*.gz"))
        if not fastq_files:
            logging.error(f"not find corresponding fastq file for {sample},please check it")
            return "NULL"
        else:
            r1 = [f for f in fastq_files if re.search(r"(_R?1)[^0-9]*\.f(ast)?q", str(f))]
            r2 = [f for f in fastq_files if re.search(r"(_R?2)[^0-9]*\.f(ast)?q", str(f))]
            if r1 and r2:
                return "PAIRED"
            else:
                return "SINGLE"            

    def classify_organism(self, df: pd.DataFrame) -> Dict[str, Dict[str, List[str]]]:
        """
        Function: Divide into different organisms, then group by single/double-ended.

        Datasheet (df) must contain the columns: sample_id, organism
        """
        groups: DefaultDict[str, DefaultDict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
        for _, row in df.iterrows():
            sample = row["sample_id"]
            organism = row["organism"]
            TYPE = self.isPairedEndSample(sample)  # 'SE' 或 'PE'
            groups[organism][TYPE].append(sample)
        return groups

    def combineFastq(self,df:pd.DataFrame):
        """
        重命名合并fastq文件
        """
        fqDir = Path(self.fqDir)
        outDir = fqDir
        samples = df['sample_id'].unique()
        for sample in samples:
            df_sample = df[df['sample_id'] == sample]
            sra_numbers = df_sample['data_id'].values
            logging.info(f"\n{sample}: {sra_numbers}")

            out_r1 = os.path.join(outDir, f"{sample}_1.fastq.gz")
            out_r2 = os.path.join(outDir, f"{sample}_2.fastq.gz")
            out_single = os.path.join(outDir, f"{sample}.fastq.gz")
            if os.path.exists(out_r1) and os.path.exists(out_r2):
                logging.info(f"[SKIP] Sample {sample} already has merged PE files.")
                continue
            if os.path.exists(out_single):
                logging.info(f"[SKIP] Sample {sample} already has merged SE file.")
                continue
                   
                # ------------------------------
            # 1个 SRA 情况 → 重命名
            # ------------------------------
            fq_all = []

            if len(sra_numbers) == 1:
                sra_number = sra_numbers[0]
                fq_files = sorted(fqDir.glob(f"{sra_number}*.gz"))
                fq_all.extend(fq_files)
                if not fq_files:
                    logging.warning(f"No fastq files found for {sample}")
                    continue
                fq_r1 = [f for f in fq_all if '_1' in os.path.basename(f)]
                fq_r2 = [f for f in fq_all if '_2' in os.path.basename(f)]
                if fq_r1 and fq_r2:
                    logging.info(f"Deteceted paired-end sample: {sample}, 2 data_id")
                    out_r1 = os.path.join(outDir, f"{sample}_1.fastq.gz")
                    out_r2 = os.path.join(outDir, f"{sample}_2.fastq.gz")
                    cmd_r1 = f"mv {fq_r1[0]} {out_r1}"
                    cmd_r2 = f"mv {fq_r2[0]} {out_r2}"
                    logging.info(f"prepare to excute command:\n{cmd_r1}\n{cmd_r2}")
                    subprocess.run(cmd_r1, shell=True, check=True)
                    subprocess.run(cmd_r2, shell=True, check=True)
                else:
                    logging.info(f"Deteceted single-end sample: {sample}, 1 data_id")
                    out = os.path.join(outDir, f"{sample}.fastq.gz")
                    cmd = f"mv {fq_all[0]} {out}"
                    logging.info(f"prepare to excute command:\n{cmd}")
                    subprocess.run(cmd, shell=True, check=True)
            # ------------------------------
            # 多个 SRA 情况 → 合并重命名
            # ------------------------------
            else:
                for sra_number in sra_numbers:
                    fq_files = sorted(fqDir.glob(f"{sra_number}*gz"))
                    fq_all.extend(fq_files)

                if not fq_all:
                    logging.warning(f"NO fastq files found for {sample}")
                    continue
                fq_r1 = [f for f in fq_all if '_1' in os.path.basename(f)]
                fq_r2 = [f for f in fq_all if '_2' in os.path.basename(f)]

                if fq_r1 and fq_r2:
                    logging.info(f"Deteceted paired-end sample: {sample}, {len(sra_numbers)} data_id")
                    fq_r1.sort()
                    fq_r2.sort()
                    out_r1 = os.path.join(outDir, f"{sample}_1.fastq.gz")
                    out_r2 = os.path.join(outDir, f"{sample}_2.fastq.gz")

                    cmd_r1 = f"cat {' '.join(fq_r1)} > {out_r1}"
                    cmd_r2 = f"cat {' '.join(fq_r2)} > {out_r2}"
                    cmd_rm1 = f"rm -f {' '.join(fq_r1)}"
                    cmd_rm2 = f"rm -f {' '.join(fq_r2)}"
                    logging.info(f"prepare excute command:\n{cmd_r1}\n{cmd_r2}\n{cmd_rm1}\n{cmd_rm2}")
                    subprocess.run(cmd_r1, shell=True, check=True)
                    subprocess.run(cmd_r2, shell=True, check=True)
                    subprocess.run(cmd_rm1,shell=True,check=True)
                    subprocess.run(cmd_rm2,shell=True,check=True)
                    logging.info(f"Merged paired files:\n  {out_r1}\n  {out_r2}")
                else:
                    logging.info(f"Detected single-end sample: {sample}, {len(sra_numbers)} data_id")
                    fq_all.sort()
                    out_single = os.path.join(outDir, f"{sample}.fastq.gz")
                    cmd = f"cat {' '.join(fq_all)} > {out_single}"
                    cmd_rm = f"rm -f {' '.join(fq_all)}"
                    logging.info(f"prepare excute command:\n{cmd}\n{cmd_rm}")
                    subprocess.run(cmd, shell=True, check=True)
                    subprocess.run(cmd_rm, shell=True, check=True)
                    logging.info(f"Merged single-end files:\n  {out_single}")
    def run(self):
        logging.info(f"load metadata from {self.meta}")
        df = self.loadMeta()
        logging.info(f"metadata was successfully load")
        try:
            logging.info(f"combine and rename fastq file to sample_id.fastq.gz")
            self.combineFastq(df)
        except Exception as e:
            logging.error(f"combine fastq file failed,error message: {e}")
        logging.info(f"classify fastq library strategy, single-end or paired-end")
        logging.info(f"To determine the species to which each sample belongs, and whether it was single-end sequencing or paired-end sequencing.")
        return self.classify_organism(df)


    


if __name__ == "__main__":
    snpMetadata = SNPMetadata("/disk5/luosg/Totipotent20251031/data/target_fq.tsv","/disk5/luosg/Totipotent20251031/data/fq")
    list1 = snpMetadata.run()
    print(list1)