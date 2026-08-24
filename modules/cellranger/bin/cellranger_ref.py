"""Build Cell Ranger reference genome from Ensembl FASTA + GENCODE GTF.

Steps (mirrors 10x Genomics GRCh38-2024-A build script):
  1. Download source FASTA/GTF (or use local paths)
  2. Modify FASTA headers: add chr prefix, handle chrM
  3. Strip Ensembl ID version suffixes from GTF
  4. Filter GTF by biotype allowlist (protein_coding, lncRNA, IG/TR genes)
  5. Remove PAR_Y genes from chrY
  6. Run cellranger mkref
"""
import argparse
import os
import re
import subprocess
import sys
import urllib.request


# 10x Genomics standard biotype allowlist
BIOTYPE_PATTERN = (
    "protein_coding|protein_coding_LoF|lncRNA|"
    "IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|"
    "IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|"
    "TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|"
    "TR_V_pseudogene|TR_J_pseudogene"
)


def download_or_copy(url_or_path, dest):
    """Download URL or copy local file to dest."""
    if url_or_path.startswith("http"):
        print(f"Downloading {url_or_path} ...")
        urllib.request.urlretrieve(url_or_path, dest)
    else:
        import shutil
        shutil.copy2(url_or_path, dest)


def modify_fasta_headers(fasta_in, fasta_out):
    """Add chr prefix to autosomes/sex chr, handle chrM.

    Input:  >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
    Output: >chr1 1
    """
    import gzip
    opener = gzip.open if fasta_in.endswith(".gz") else open
    with opener(fasta_in, "rt") as fin, open(fasta_out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                # Replace metadata after space with contig name
                parts = line[1:].split(None, 1)
                name = parts[0]
                chr_name = name
                if re.match(r"^[0-9]+$", name) or name in ("X", "Y"):
                    chr_name = f"chr{name}"
                elif name == "MT":
                    chr_name = "chrM"
                line = f">{chr_name} {name}\n"
            fout.write(line)


def modify_gtf_ids(gtf_in, gtf_out):
    """Strip version suffixes from Ensembl gene/transcript/exon IDs.

    Input:  gene_id "ENSG00000223972.5";
    Output: gene_id "ENSG00000223972"; gene_version "5";
    """
    id_pattern = re.compile(
        r'(ENS(?:MUS)?[GTE]\d+)\.(\d+)'
    )
    import gzip
    opener = gzip.open if gtf_in.endswith(".gz") else open
    with opener(gtf_in, "rt") as fin, open(gtf_out, "w") as fout:
        for line in fin:
            for id_type in ("gene_id", "transcript_id", "exon_id"):
                line = re.sub(
                    rf'{id_type} "(ENS(?:MUS)?[GTE]\d+)\.(\d+)";',
                    rf'{id_type} "\1"; {id_type.replace("_id", "_version")} "\2";',
                    line,
                )
            fout.write(line)


def filter_gtf_by_biotype(gtf_in, gtf_out, biotype_re):
    """Filter GTF: keep only allowed biotypes, exclude readthrough, remove PAR_Y."""
    gene_pattern = re.compile(rf'gene_type "({BIOTYPE_PATTERN})"')
    tx_pattern = re.compile(rf'transcript_type "({BIOTYPE_PATTERN})"')
    readthrough_pattern = re.compile(r'tag "readthrough_transcript"')
    gene_id_pattern = re.compile(r'(gene_id "[^"]+")')

    # Step 1: collect gene IDs from transcripts passing filters
    gene_ids = set()
    import gzip
    opener = gzip.open if gtf_in.endswith(".gz") else open
    with opener(gtf_in, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            if not gene_pattern.search(attrs):
                continue
            if not tx_pattern.search(attrs):
                continue
            if readthrough_pattern.search(attrs):
                continue
            m = gene_id_pattern.search(attrs)
            if m:
                gene_ids.add(m.group(1))

    print(f"  {len(gene_ids)} genes pass biotype filter")

    # Step 2: write header
    with open(gtf_in, "r") as fin, open(gtf_out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
            else:
                break

    # Step 3: filter to allowlisted genes, remove PAR_Y
    with open(gtf_in, "r") as fin, open(gtf_out, "a") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            # Check gene is in allowlist
            if not any(gid in fields[8] for gid in gene_ids):
                continue
            # Remove PAR_Y: exclude chrY entries in PAR range (except ENSG00000290840)
            if fields[0] == "chrY":
                start = int(fields[3])
                if start < 2752083 or start >= 56887903:
                    continue
                if "ENSG00000290840" in fields[8]:
                    continue
            fout.write(line)


def run_cellranger_mkref(cellranger, genome_name, version, fasta, genes, nthreads):
    """Run cellranger mkref."""
    cmd = [
        cellranger, "mkref",
        f"--ref-version={version}",
        f"--genome={genome_name}",
        f"--fasta={fasta}",
        f"--genes={genes}",
        f"--nthreads={nthreads}",
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.check_call(cmd)


def main():
    parser = argparse.ArgumentParser(
        description="Build Cell Ranger reference from FASTA + GTF"
    )
    parser.add_argument("--fasta", required=True, help="Ensembl FASTA URL or local path")
    parser.add_argument("--gtf", required=True, help="GENCODE GTF URL or local path")
    parser.add_argument("--output", required=True, help="Output reference directory")
    parser.add_argument("--genome", default="GRCh38", help="Genome name")
    parser.add_argument("--version", default="2024-A", help="Reference version")
    parser.add_argument("--cellranger", default="cellranger", help="Cell Ranger binary path")
    parser.add_argument("--nthreads", type=int, default=16, help="Threads for mkref")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    build_dir = os.path.join(args.output, "build")
    source_dir = os.path.join(build_dir, "source")
    os.makedirs(source_dir, exist_ok=True)

    # Derive filenames
    fasta_basename = os.path.basename(args.fasta).replace(".gz", "")
    gtf_basename = os.path.basename(args.gtf).replace(".gz", "")
    fasta_src = os.path.join(source_dir, fasta_basename)
    gtf_src = os.path.join(source_dir, gtf_basename)
    fasta_modified = os.path.join(build_dir, f"{fasta_basename}.modified")
    gtf_modified = os.path.join(build_dir, f"{gtf_basename}.modified")
    gtf_filtered = os.path.join(build_dir, f"{gtf_basename}.filtered")

    # Step 1: Download
    if not os.path.exists(fasta_src):
        download_or_copy(args.fasta, fasta_src)
    else:
        print(f"  FASTA already exists: {fasta_src}")

    if not os.path.exists(gtf_src):
        download_or_copy(args.gtf, gtf_src)
    else:
        print(f"  GTF already exists: {gtf_src}")

    # Step 2: Modify FASTA headers
    print("Modifying FASTA headers ...")
    modify_fasta_headers(fasta_src, fasta_modified)

    # Step 3: Modify GTF IDs
    print("Modifying GTF IDs ...")
    modify_gtf_ids(gtf_src, gtf_modified)

    # Step 4: Filter GTF
    print("Filtering GTF by biotype ...")
    filter_gtf_by_biotype(gtf_modified, gtf_filtered, BIOTYPE_PATTERN)

    # Step 5: Build reference
    print("Building Cell Ranger reference ...")
    run_cellranger_mkref(
        args.cellranger, args.genome, args.version,
        fasta_modified, gtf_filtered, args.nthreads,
    )

    print(f"Done. Reference built at {args.output}")


if __name__ == "__main__":
    main()
