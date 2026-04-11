import ftplib
import logging
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.LogUtil import setup_logger
from common.SepUtil import detect_delimiter
from common.MatchUtil import run_accession_match
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, List
import logging
logger = logging.getLogger(__name__)

def download_files_from_ftp(server, username=None, password=None, remote_dir="/", local_dir="./downloads"):
    """
    Download all files from a specified directory on an FTP server to a local directory.

    :param server: FTP server address
    :param username: FTP username (optional, default is None for anonymous login)
    :param password: FTP password (optional, default is None for anonymous login)
    :param remote_dir: Directory on the FTP server to download files from
    :param local_dir: Local directory to save the downloaded files
    """
    try:
        # Check if remote_dir contains the server address
        if remote_dir.startswith("ftp://"):
            remote_dir = remote_dir.replace(f"ftp://{server}", "")

        logger.info(f"Connecting to FTP server: {server}")
        ftp = ftplib.FTP(server)
        if username and password:
            ftp.login(user=username, passwd=password)
            logger.info("Logged in with provided credentials.")
        else:
            ftp.login()  # Anonymous login
            logger.info("Logged in anonymously.")

        # Change to the specified remote directory
        logger.info(f"Changing to remote directory: {remote_dir}")
        ftp.cwd(remote_dir)

        # Ensure the local directory exists
        if not os.path.exists(local_dir):
            os.makedirs(local_dir)
            logger.info(f"Created local directory: {local_dir}")

        # List files in the remote directory
        files = ftp.nlst()
        logger.info(f"Files found in remote directory: {files}")

        for file in files:
            local_file_path = os.path.join(local_dir, file)
            with open(local_file_path, 'wb') as local_file:
                ftp.retrbinary(f'RETR {file}', local_file.write)
                logger.info(f"Downloaded: {file} to {local_file_path}")

        # Close the FTP connection
        ftp.quit()
        logger.info("All files downloaded successfully.")

    except ftplib.all_errors as e:
        logger.error(f"FTP error: {e}")

def download_file_from_ftp(server, username=None, password=None, remote_file_path="", local_file_path=""):
    """
    Download a specific file from an FTP server.

    :param server: FTP server address
    :param username: FTP username (optional, default is None for anonymous login)
    :param password: FTP password (optional, default is None for anonymous login)
    :param remote_file_path: Path to the file on the FTP server
    :param local_file_path: Path to save the file locally
    """
    try:
        # Check if remote_file_path contains the server address
        if remote_file_path.startswith("ftp://"):
            remote_file_path = remote_file_path.replace(f"ftp://{server}", "")

        logger.info(f"Connecting to FTP server: {server}")
        ftp = ftplib.FTP(server)
        if username and password:
            ftp.login(user=username, passwd=password)
            logger.info("Logged in with provided credentials.")
        else:
            ftp.login()  # Anonymous login
            logger.info("Logged in anonymously.")

        # Ensure the local directory exists
        local_dir = os.path.dirname(local_file_path)
        if local_dir and not os.path.exists(local_dir):
            os.makedirs(local_dir)
            logger.info(f"Created local directory: {local_dir}")

        # Download the specific file
        with open(local_file_path, 'wb') as local_file:
            ftp.retrbinary(f'RETR {remote_file_path}', local_file.write)
            logger.info(f"Downloaded: {remote_file_path} to {local_file_path}")

        # Close the FTP connection
        ftp.quit()
        logger.info("File downloaded successfully.")

    except ftplib.all_errors as e:
        logger.error(f"FTP error: {e}")

def parser_metadata(
    meta_experiment_file:str,
    ftp_ref_file:str,
    output_dir:str,
    condition_col:str = "experiment_title",
    condition_value:str = "Cell line co-culture (mouse and human ips)",
    run_accession_col:str = "run_accession",
):
    """
    Function: Parse metadata from downloaded files(CNGBdb).
    Parameters:
    - meta_experiment_file: Path to the metadata file containing experiment information.
    - ftp_ref_file: Path to the reference file containing FTP URLs. url should be in the second column.
    - condition_col: Column name in the metadata file to filter by condition (default: "Experiment Title").
    - condition_value: Value in the condition column to filter the metadata (default: "Cell line co-culture (mouse and human ips)").
    - project_accession_col: Column name for project accession in the metadata file (default: "project_accession").
    - sample_accession_col: Column name for sample accession in the metadata file (default: "sample_accession").
    - experiment_accession_col: Column name for experiment accession in the metadata file (default: "experiment_accession").
    - run_accession_col: Column name for run accession in the metadata file (default: "run_accession").
    Returns:
    - A filtered DataFrame containing metadata records that match the specified condition.
    - A list of FTP URLs extracted from the reference file.
    """
    meta_sep = detect_delimiter(meta_experiment_file)
    df = pd.read_csv(meta_experiment_file, sep=meta_sep)
    df_filtered = df[df[condition_col] == condition_value]
    logger.info(f"Filtered metadata: {len(df_filtered)} records match the condition '{condition_value}' in column '{condition_col}'.")
    ftp_urls = []

    with open(ftp_ref_file, 'r', encoding='utf-8') as f:
        for line in f:
            ftp_urls.append(line.strip().split()[1])
    run_accessions = df_filtered[run_accession_col].tolist()
    matched_urls, unmatched_accessions = run_accession_match(ftp_urls, run_accessions)
    
    with open(os.path.join(output_dir, "unmatched_accessions.txt"), 'w', encoding='utf-8') as f:
        for acc in unmatched_accessions:
            f.write("\t".join(acc) + "\n")
    
    with open(os.path.join(output_dir, "matched_urls.txt"), 'w', encoding='utf-8') as f:
        for url in matched_urls:
            f.write(url + "\n")

    return matched_urls

def main():
    experiment_metadata_file = "/data/pub/zhousha/20260411_RNAseq/data/meta/metadata_CNP0003135_experiment.tsv"
    ftp_ref_file = "/data/pub/zhousha/20260411_RNAseq/data/meta/data_download_links_CNP0003135_ftp.txt"
    output_dir = "/data/pub/zhousha/20260411_RNAseq/data/fastq"
    log_file = "/data/pub/zhousha/20260411_RNAseq/log/ftp_download.log"
    logger = setup_logger("", log_file=log_file)
    ftp_urls = parser_metadata(
        meta_experiment_file=experiment_metadata_file,
        ftp_ref_file=ftp_ref_file,
        output_dir=output_dir
    )
    logger.info(f"Total matched FTP URLs: {len(ftp_urls)}")
    # Parallel download for FTP URLs
    def download_url(url):
        if url.startswith("ftp://"):
            server = url.split("/")[2]
            remote_file_path = "/" + "/".join(url.split("/")[3:])
            local_file_path = os.path.join(output_dir, os.path.basename(remote_file_path))
            logger.info(f"Server: {server}, Remote file path: {remote_file_path}, Local file path: {local_file_path}")
            download_file_from_ftp(server, remote_file_path=remote_file_path, local_file_path=local_file_path)

    with ThreadPoolExecutor(max_workers=4) as executor:
        executor.map(download_url, ftp_urls)

    logger.info("All downloads completed.")

if __name__ == "__main__":
    main()