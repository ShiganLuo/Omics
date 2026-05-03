from typing import List
import argparse
import logging
import glob

logging.basicConfig(level=logging.INFO)
logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

def ucsc_track_format(
    bigwig_files:List[str],
    output:str
):
    """
    Convert a list of bigwig files to UCSC track format.

    Args:
        bigwig_files (List[str]): A list of bigwig file paths.

    Returns:
        str: A string in UCSC track format.
    """
    track_lines = []
    for bigwig_file in bigwig_files:
        track_name = bigwig_file.split('/')[-1].replace('.bigwig', '')
        track_line = f'track type=bigWig name="{track_name}" description="{track_name}" bigDataUrl={bigwig_file}'
        track_lines.append(track_line)
    with open(output, 'w') as f:
        f.write('\n'.join(track_lines))
    return track_lines

def test_ucsc_track_format():
    bigwig_files = glob.glob('/data/pub/zhousha/20260422_ClIPseq/output/CLIP/igv/*.bigwig')
    output = "/data/pub/zhousha/20260422_ClIPseq/output/CLIP/igv/ucsc_tracks.txt"
    ucsc_track_format(bigwig_files, output)

def main():
    parser = argparse.ArgumentParser(description='Convert bigwig files to UCSC track format.')
    parser.add_argument('-i', '--bigwig_files', required=True, nargs='+', help='A list of bigwig file paths.')
    parser.add_argument('-o', '--output', required=True, help='Output file to save the UCSC track format.')
    args = parser.parse_args()
    logging.info(f"Converting {len(args.bigwig_files)} bigwig files to UCSC track format.")
    track_format = ucsc_track_format(args.bigwig_files, args.output)
    logging.info("Conversion complete.")

if __name__ == "__main__":
    # main()
    test_ucsc_track_format()