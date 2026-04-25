import argparse
import pandas as pd
from typing import List

def combine_files_by_line(input_files: List[str], output_file: str) -> None:
    """Combine multiple files by line and save to a new file."""
    combined_df = pd.DataFrame()
    df_list = []
    for file in input_files:
        df = pd.read_csv(file,sep="\t")
        df_list.append(df)
    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df = combined_df.sort_values()
    combined_df.to_csv(output_file, index=False, sep="\t")

def main():
	parser = argparse.ArgumentParser(description="Combine files by line using pandas.")
	parser.add_argument('-i', '--input', nargs='+', required=True, help='Input files')
	parser.add_argument('-o', '--output', required=True, help='Output file')
	args = parser.parse_args()
	combine_files_by_line(args.input, args.output)

if __name__ == '__main__':
	main()
