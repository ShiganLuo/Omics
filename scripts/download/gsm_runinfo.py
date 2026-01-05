import os
import re
import csv
from time import sleep
from typing import List, Dict
import subprocess
import logging
import random
import pandas as pd
import xml.etree.ElementTree as ET
from typing import Callable, List
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	datefmt='%Y-%m-%d %H:%M:%S'

)
def extract_geo_info(filepath, output_tsv):
    """
    Extract all entries (Title, Organism, FTP, Accession) from a GEO text file, even if the order or line breaks are not fixed, and save it as TSV.

    Args:

    Returns:
        A list of dictionaries containing the extraction results for each GSM record.

        {
            "Accession": accession,
            "Title": title,
            "FTP_Download": ftp,
            "Organism": organism
        }
    """
    results = []

    if not os.path.exists(filepath):
        logging.error(f"file not found:{filepath}")
        return results

    with open(filepath, 'r', encoding='utf-8') as f:
        text = f.read()

    # Split each entry into chunks
    entry_blocks = re.split(r'^\s*\d+\.', text, flags=re.MULTILINE)
    for block in entry_blocks:

        title_match = re.match(r'(.+)', block)
        title = title_match.group(1).strip() if title_match else "N/A"

        organism_match = re.search(r'Organism:\s*([^\n]+)', block)
        organism = organism_match.group(1).strip() if organism_match else "N/A"

        ftp_match = re.search(r'FTP download:\s*.*(ftp://[^\s]+)', block)
        ftp = ftp_match.group(1).strip() if ftp_match else "N/A"

        accession_match = re.search(r'Series\s+Accession:\s*(\S+)', block)
        accession = accession_match.group(1).split('\t')[0] if accession_match else "N/A"
        if accession == "N/A":
            continue
        results.append({
            "Accession": accession,
            "Title": title,
            "FTP_Download": ftp,
            "Organism": organism
        })

    if results:
        keys = ["Accession", "Title", "FTP_Download", "Organism"]
        with open(output_tsv, 'w', newline='', encoding='utf-8-sig') as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=keys, delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
        logging.info(f"已保存 {len(results)} 条记录到 {output_tsv}")

    return results

def parse_geo_summary_json(
        filepath: str,
        output_tsv: str,
        write_mode: str = 'w'
) -> List[Dict[str, str]]:
    """
    Parse the GEO summary JSON data and extract key information from the block where Accession is GSM.
    
    Args:
        filepath: the path of outputfile by `esearch -db gds -query GSE308214 |  efetch -format json` command
        
    Returns:
        A list of dictionaries containing the extraction results for each GSM record.
        {
            "Accession": GSM……,
            "Series": GSE……,
            "FTP download": ftp,
            "Organism": organism
            "Source name": source_name
        }
    """
    results = []
    if not os.path.exists(filepath):
        logging.error(f"file not found {filepath}")
        return results

    with open(filepath, 'r', encoding='utf-8') as f:
        text = f.read()
        
    patterns = {
        "Accession": r'Accession:\s*(GSM\d+)',
        "FTP download": r'FTP download:.*?ftp://(.*?)(?:\s|$)',
        "Series": r'Series:.*?(GSE\d+)',
        "Organism": r'Organism:\s*(.*?)\n',
        "Source name": r'Source name:\s*(.*?)\n',
    }
    entry_blocks = re.split(r'^\s*\d+\.', text, flags=re.MULTILINE)
    for block in entry_blocks:
        if "Accession: GSM" in block:
            data = {}
            for field, pattern in patterns.items():
                match = re.search(pattern, block, re.DOTALL | re.IGNORECASE)
                if match:
                    if field == "FTP download":
                        ftp_match = re.search(r'FTP download:.*?ftp://(\S+)', block, re.DOTALL | re.IGNORECASE)
                        data[field] = 'ftp://' + ftp_match.group(1).strip() if ftp_match else 'N/A'
                    elif field == "Series":
                        data[field] = match.group(1).strip()
                    else:
                        data[field] = match.group(1).strip()
                else:
                    data[field] = 'N/A'

            if data.get("Accession", "N/A").startswith("GSM"):
                final_accession_match = re.search(r'Accession:\s*(GSM\d+)', block, re.IGNORECASE)
                if final_accession_match:
                    data['Accession'] = final_accession_match.group(1).strip()

                results.append(data)
    if results:
        file_exists = os.path.exists(output_tsv) and os.path.getsize(output_tsv) > 0
        keys = ["Accession", "Series", "FTP download", "Organism", "Source name"]
        with open(output_tsv, write_mode, newline='', encoding='utf-8-sig') as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=keys, delimiter='\t')
            if not file_exists:
                writer.writeheader()
            writer.writerows(results)
        logging.info(f"已保存 {len(results)} 条记录到 {output_tsv}")
    return results

def parse_geo_summary_xml(
        file_path:str,
        output_tsv:str,
        write_mode:str='w'
) -> List[Dict[str, str]]:
    results = []
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
    except ET.ParseError as e:
        logging.error(f"XML解析错误:{e}, {file_path}")
        return  results
    except FileNotFoundError:
        logging.error(f"文件未找到：{file_path}")
        return results
    
    for doc_summary in root.findall('DocumentSummary'):
        accession = doc_summary.find('Accession').text if doc_summary.find('Accession').text is not None else "N/A"
        if "GSE" in accession:
            gse_id = accession
            title = doc_summary.find('title').text if doc_summary.find('title') is not None and doc_summary.find('title').text is not None else "N/A"
            gds_type = doc_summary.find('gdsType').text if doc_summary.find('gdsType') is not None and doc_summary.find('gdsType').text is not None else "N/A"
            publish_date = doc_summary.find('PDAT').text if doc_summary.find('PDAT') is not None and doc_summary.find('PDAT').text is not None else "N/A"
            n_samples = doc_summary.find('n_samples').text if doc_summary.find('n_samples') is not None and doc_summary.find('n_samples').text is not None else "N/A"
            FTPLink = doc_summary.find('FTPLink').text if doc_summary.find('FTPLink') is not None and doc_summary.find('FTPLink').text is not None else "N/A"
            BioProject = doc_summary.find('BioProject').text if doc_summary.find('BioProject') is not None and doc_summary.find('BioProject').text is not None else "N/A"
            taxon = doc_summary.find('taxon').text if doc_summary.find('taxon') is not None and doc_summary.find('taxon').text is not None else "N/A"
            summary = doc_summary.find('summary').text if doc_summary.find('summary') is not None and doc_summary.find('summary').text is not None else "N/A"
            
            for sample_element in doc_summary.findall('./Samples/Sample'):
                gsm_idFromGSE = sample_element.find('Accession').text if sample_element.find('Accession') is not None and sample_element.find('Accession').text is not None else "N/A"
                gsm_titleFromGSE = sample_element.find('Title').text if sample_element.find('Title') is not None and sample_element.find('Title').text is not None else "N/A"
                results.append({
                    "GSMId":gsm_idFromGSE,
                    "GSMDescription": gsm_titleFromGSE,
                    "GSEId": gse_id,
                    "GSEDescription":title,
                    "GSEDetailedDescription":summary,
                    "Organism":taxon,
                    "GDSType": gds_type,
                    "PublishDate":publish_date,
                    "SamplesNum":n_samples,
                    "BioProject":BioProject,
                    "FTPLink":FTPLink
                })
        elif "GDS" in accession:
            gds_id = accession
            continue
        elif "GPL" in accession:
            gpl_id = accession
            continue
        elif "GSM" in accession:
            gsm_id = accession
            continue
        else:
            logging.warning(f"出现未预料的accession,文件: {file_path}")
            continue
    if results:
        file_exists = os.path.exists(output_tsv) and os.path.getsize(output_tsv) > 0
        keys = ["GSMId", "GSMDescription", "GSEId", "GSEDescription","GSEDetailedDescription","Organism", "GDSType","PublishDate","SamplesNum","BioProject","FTPLink"]
        with open(output_tsv, write_mode, newline='', encoding='utf-8-sig') as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=keys, delimiter='\t')
            if not file_exists:
                writer.writeheader()
            writer.writerows(results)
        logging.info(f"已保存 {len(results)} 条记录到 {output_tsv}")
    return results
            
    

def batch_extract_info(folder_path:str,
                       output_tsv:str,
                       suffix:str=".json",
                       callback: Callable[[str,str,str],List] = parse_geo_summary_json):
    """
    Batch process all .json files in a folder.
    """
    all_results = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith(suffix):
            filepath = os.path.join(folder_path, filename)
            logging.info(f"正在处理文件: {filepath}")
            results = callback(filepath,output_tsv,write_mode = 'a')
            all_results.extend(results)
            
    return all_results

def runEdirectCommand(id:str,outfile:str,database:str="gds",format:str="json"):
    """
    excute command: `esearch -db [database] -query [id] |  efetch -format [format] > outfile` for GSEid
    format: json,docsum,runinfo
    database: gds,sra
    """
    with open(outfile, "w") as f_out:
        p1 = subprocess.Popen(["esearch", "-db", database, "-query", id], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        p2 = subprocess.Popen(["efetch", "-format", format], stdin=p1.stdout, stdout=f_out,stderr=subprocess.PIPE)
        p1.stdout.close()  # 关闭父进程对p1输出的引用
        p2.communicate()   # 等待p2完成
    if os.path.getsize(outfile) == 0:
        logging.error(f"the {outfile} has no content")
    n = random.randint(0, 2)
    sleep(n)
    # subprocess.run(Command)



def removeDuplicate(filepath:str,output_tsv:str):
    """
    remove duplicate line according to the first column
    """
    df = pd.read_csv(filepath,sep="\t",header=None,keep_default_na=False)
    df = df.drop_duplicates(0,keep="first")
    df.to_csv(output_tsv,sep="\t",index=False,header=False)

def runinfo_formated(folder_path:str,outfile:str,suffix:str,write_mode:str="a"):
    """
    extract needed information from runinfo
    """
    
    for filename in os.listdir(folder_path):
        if filename.endswith(suffix):
            filepath = os.path.join(folder_path, filename)
            logging.info(f"正在处理文件: {filepath}")
        
            isInfileExists = os.path.exists(filepath) and os.path.getsize(filepath) > 0
            if not isInfileExists:
                logging.warning(f"{filepath} has not content, jump over")
                continue
            try:
                df = pd.read_csv(filepath,sep=",")
            except Exception as e:
                logging.error(f"{filepath} read filed, error: {e}")
                continue
            df = df[["Run","SampleName","BioProject","ScientificName","Sex","Disease","Tumor","CenterName","LibraryStrategy","LibraryLayout","avgLength"]]
            file_exists = os.path.exists(outfile) and os.path.getsize(outfile) > 0
            with open(outfile,write_mode) as f:
                df.to_csv(f,sep="\t",index=False,header=not file_exists)



if __name__ == '__main__':
    # 调用函数
    # # extrat formated information from gds result
    file_path = '/home/luosg/Data/genomeStability/data/gds_result.txt' # 将此处的 'data.txt' 替换为你的实际文件名
    output_tsv = "/home/luosg/Data/genomeStability/data/GEO_related.tsv"
    # results = extract_geo_info(file_path,output_tsv)

    #### get GSE json information
    # outdir = "/home/luosg/Data/genomeStability/data/json"
    # for GSE_dict in results:
    #     GSEId = GSE_dict["Accession"]
    #     outfile = f"{outdir}/{GSEId}.json"
    #     logging.info(f"process {GSEId}, output it to {outfile}")
    #     runEdirectCommand(GSEId,outfile,"json")

    #### get GSE xml information
    # outdir = "/home/luosg/Data/genomeStability/data/html"
    # for GSE_dict in results:
    #     GSEId = GSE_dict["Accession"]
    #     outfile = f"{outdir}/{GSEId}.html"
    #     logging.info(f"process {GSEId}, output it to {outfile}")
    #     runEdirectCommand(GSEId,outfile,"docsum")

    #### get GSE xml information (wrong download in last step)
    # outdir = "/home/luosg/Data/genomeStability/data/html"
    # df = pd.read_csv("/home/luosg/Data/genomeStability/log/xml_wrong.log",header=None)
    # for GSEId in df[0]:
    #     outfile = f"{outdir}/{GSEId}.html"
    #     logging.info(f"process {GSEId}, output it to {outfile}")
    #     runEdirectCommand(GSEId,outfile,"docsum")

    
    #### extract formated information from GSE json
    # folder_path = "/home/luosg/Data/genomeStability/data/json"
    # output_tsv = "/home/luosg/Data/genomeStability/data/GSM_FromJson.tsv"
    # batch_extract_info(folder_path,output_tsv,'.json',parse_geo_summary_json)

    #### remove GSM_FromJson duplicates
    # removeDuplicate(output_tsv,"/home/luosg/Data/genomeStability/data/GSM_FromJson_dedup.tsv")

    ###extract formated information from GSE xml
    # folder_path = "/home/luosg/Data/genomeStability/data/html"
    output_tsv = "/home/luosg/Data/genomeStability/data/GSE_FromXml.tsv"
    # batch_extract_info(folder_path,output_tsv,'.html',parse_geo_summary_xml)
    removeDuplicate(output_tsv,"/home/luosg/Data/genomeStability/data/GSE_FromXml_dedup.tsv")

    # get sra runinfo according to PRJId
    # df = pd.read_csv("/home/luosg/Data/genomeStability/data/GSE_FromXml_dedup_simplify.tsv",sep="\t",header=None)
    # PRJIds = df[4].unique()
    # outdir = "/home/luosg/Data/genomeStability/data/sra"
    # for PRJid in PRJIds:
    #     outfile = f"{outdir}/{PRJid}_runinfo.csv"
    #     logging.info(f"process {PRJid}, output it to {outfile}")
    #     runEdirectCommand(PRJid,outfile,"sra","runinfo")
    #     # runEdirectCommand(PRJid,outfile,"sra","runinfo")

    # combine sra runinfo information
    # folder_path = "/home/luosg/Data/genomeStability/data/sra"
    # output_tsv = "/home/luosg/Data/genomeStability/data/SRA_runinfo.tsv"
    # runinfo_formated(folder_path,output_tsv,"_runinfo.csv","a")
    # removeDuplicate(output_tsv,"/home/luosg/Data/genomeStability/data/SRA_runinfo_dedup.tsv")




    

    

