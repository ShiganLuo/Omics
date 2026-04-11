from typing import List,Tuple
import logging
logger = logging.getLogger(__name__)
def run_accession_match(
    ftp_urls: List[str],
    run_accessions: List[str]
)-> Tuple[List[str], List[str]]:
    """
    Match metadata accessions with FTP URLs.

    Parameters:
    - ftp_urls: List of FTP URLs to search through.
    - run_accessions: List of run accessions to match.

    Returns:
    - A list of matched FTP URLs that contain all specified accessions.
    - A list of unmatched accessions for logging or further analysis.
    """
    matched_urls = []
    unmatched_accessions = []
    for run_acc in run_accessions:
        matched = [p for p in ftp_urls if run_acc in p]
        if not matched:
            logger.warning(f"No match found for run accession '{run_acc}' in FTP URLs.")
            unmatched_accessions.append(run_acc)
        if len(matched) > 1:
            logger.warning(f"Multiple matches found for run accession '{run_acc}': {matched}")
        matched_urls.extend(matched)
    return matched_urls, unmatched_accessions
