import sys
from snakemake import config

def main():
    output_path = sys.argv[1] if len(sys.argv) > 1 else 'group.tsv'
    control_samples = config.get('control_samples', [])
    treatment_samples = config.get('treatment_samples', [])
    control_group_name = config.get('control_group_name', 'control')
    treatment_group_name = config.get('treatment_group_name', 'treatment')
    with open(output_path, 'w') as f:
        f.write('sample\tgroup\n')
        for s in control_samples:
            f.write(f'{s}\t{control_group_name}\n')
        for s in treatment_samples:
            f.write(f'{s}\t{treatment_group_name}\n')

if __name__ == '__main__':
	main()
