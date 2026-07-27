#!/usr/bin/env python3

import argparse


def _split_samples(value: str) -> list[str]:
    if not value:
        return []
    return [item for item in value.split(",") if item]


def main() -> None:
    parser = argparse.ArgumentParser(description="Write group.tsv for DESeq2 analysis")
    parser.add_argument("-o", "--output", default="group.tsv", help="Output group TSV path")
    parser.add_argument("-c", "--control-samples", default="", help="Comma-separated control sample IDs")
    parser.add_argument("-t", "--treatment-samples", default="", help="Comma-separated treatment sample IDs")
    parser.add_argument("-p", "--control-group-name", default="control", help="Control group label")
    parser.add_argument("-e", "--treatment-group-name", default="treatment", help="Treatment group label")
    args = parser.parse_args()

    control_samples = _split_samples(args.control_samples)
    treatment_samples = _split_samples(args.treatment_samples)

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write("sample\tgroup\n")
        for sample_id in control_samples:
            handle.write(f"{sample_id}\t{args.control_group_name}\n")
        for sample_id in treatment_samples:
            handle.write(f"{sample_id}\t{args.treatment_group_name}\n")


if __name__ == "__main__":
    main()
