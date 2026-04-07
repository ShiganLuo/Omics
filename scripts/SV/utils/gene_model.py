#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module: gene_model_plot_vcf
Features:
- Function-based, reusable, modular
- Draw gene model with exons/introns + mutations/SVs
- Supports GTF/GFF3 input
- Mutations input: CSV or VCF
"""

import os
from typing import List, Tuple, Dict, Optional, Union
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import gffutils
import argparse
import vcf  # pip install PyVCF

# -----------------------------
# 1. Gene structure parsing
# -----------------------------
def create_db(gtf_file: str) -> gffutils.FeatureDB:
    db_file = gtf_file + ".db"
    if not os.path.exists(db_file):
        gffutils.create_db(
            gtf_file,
            dbfn=db_file,
            force=True,
            keep_order=True,
            merge_strategy='merge',
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
    return gffutils.FeatureDB(db_file)

def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping or adjacent intervals."""
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged

<<<<<<< HEAD
def get_gene_structure_by_transcript(db: gffutils.FeatureDB, gene_name: str) -> Tuple[Dict[str, Dict[str, List[Tuple[int, int]]]], gffutils.Feature]:
=======
def get_gene_structure_by_transcript(
        db: gffutils.FeatureDB, 
        gene_name: str
    ) -> Tuple[Dict[str, Dict[str, List[Tuple[int, int]]]], gffutils.Feature]:
>>>>>>> 0ea0c995979175199379f89ae752af4c876178c6
    """
    Return gene structure grouped by transcript (all isoforms).
    Output format:
    {
        transcript_id: {
            'exons': [(start, end), ...],
            'UTRs': [(start, end), ...]
        },
        ...
    }
    """
    gene: Optional[gffutils.Feature] = None
    for g in db.features_of_type('gene'):
        if g.attributes.get('gene_name', [g.id])[0] == gene_name or g.id == gene_name:
            gene = g
            break
    if gene is None:
        raise ValueError(f"Gene {gene_name} not found in GTF")

    transcripts: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}
    for tx in db.children(gene, featuretype=('mRNA','transcript'), order_by='start'):
        tx_id: str = tx.id
        exons: List[Tuple[int,int]] = [(e.start, e.end) for e in db.children(tx, featuretype='exon', order_by='start')]
        exons = merge_intervals(exons)
        utrs: List[Tuple[int,int]] = [(u.start, u.end) for u in db.children(tx, featuretype='UTR', order_by='start')]
        utrs = merge_intervals(utrs)
        transcripts[tx_id] = {'exons': exons, 'UTRs': utrs}

    if not transcripts:
        raise ValueError(f"No transcripts found for gene {gene_name}. Check GTF feature types (mRNA/transcript).")

    return transcripts, gene

# -----------------------------
# 2. Load mutations (CSV or VCF)
# -----------------------------
Mutation = Dict[str, Union[int, str]]

def load_mutations(mutation_file: Optional[str]) -> List[Mutation]:
    if not mutation_file:
        return []
    
    if mutation_file.lower().endswith(".vcf"):
        vcf_reader = vcf.Reader(filename=mutation_file)
        mutations: List[Mutation] = []
        for record in vcf_reader:
            if record.is_snp:
                mut_type = "SNV"
            elif record.is_indel:
                mut_type = "INDEL"
            elif record.is_sv:
                mut_type = record.INFO.get('SVTYPE', 'SV')
            else:
                mut_type = "MUT"
            mutations.append({
                "pos": record.POS,
                "type": mut_type,
                "label": str(record.ID) if record.ID else ""
            })
        return mutations
    
    # fallback to CSV
    df = pd.read_csv(mutation_file)
    required_cols = ["pos","type"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' missing in mutations CSV")
    mutations = []
    for _, row in df.iterrows():
        mutations.append({
            "pos": row["pos"],
            "type": row["type"],
            "label": row.get("label", "")
        })
    return mutations

# -----------------------------
# 3. Plotting function (all transcripts)
# -----------------------------
def plot_gene_model_all_transcripts(
    gene: gffutils.Feature,
    transcripts: Dict[str, Dict[str, List[Tuple[int,int]]]],
    mutations: Optional[List[Mutation]] = None,
    figsize: Tuple[float,float]=(10,2),
    output: str="gene_model.png",
    show_mutation_label: bool=False
) -> None:
    mutations = mutations or []
    fig, ax = plt.subplots(figsize=figsize)
    n_transcripts = len(transcripts)
    y_gap = 1.0 / max(n_transcripts, 1)
    
    type_colors = {"DEL":"#0AF0A7", "DUP":"#F00A0AFA", "INS":"#0A9CF0",
                   "INV":"#9FF00A", "BND":"#E4F00A", "TRA": "#F0660A"}

    for i, (tx_id, tx_data) in enumerate(transcripts.items()):
        y_base = i * y_gap
        exons = tx_data['exons']
        utrs = tx_data['UTRs']

        # Draw introns as lines
        for j in range(len(exons)-1):
            ax.plot([exons[j][1], exons[j+1][0]], [y_base, y_base], color='black', lw=1, zorder=0)

        # Exons
        for start, end in exons:
<<<<<<< HEAD
            ax.add_patch(patches.Rectangle((start, y_base-0.1), end-start, 0.2, color='red', ec='red', zorder=2))
        
        # UTRs
        for start, end in utrs:
            ax.add_patch(patches.Rectangle((start, y_base-0.05), end-start, 0.1, color='green', ec='green', zorder=3))
=======
            ax.add_patch(patches.Rectangle((start, y_base-0.1), end-start, 0.2, color='skyblue', ec='red', zorder=2))
        
        # UTRs
        for start, end in utrs:
            ax.add_patch(patches.Rectangle((start, y_base-0.05), end-start, 0.1, color='lightgreen', ec='green', zorder=3))
>>>>>>> 0ea0c995979175199379f89ae752af4c876178c6
        
        # Transcript label
        ax.text(gene.start-50, y_base, tx_id, fontsize=8, ha='right', va='center')

    # Strand arrow
    arrow_gap = (gene.end-gene.start)/20
    if gene.strand == '+':
        ax.arrow(gene.start, -0.1, gene.end-gene.start, 0, head_width=0.05, head_length=arrow_gap,
                 length_includes_head=True, color='black')
    elif gene.strand == '-':
        ax.arrow(gene.end, -0.1, gene.start-gene.end, 0, head_width=0.05, head_length=arrow_gap,
                 length_includes_head=True, color='black')

    # Mutations/SVs
    for mut in mutations:
        color = type_colors.get(mut["type"], "gray")
        ax.scatter(mut["pos"], 1.0, color=color, s=50, zorder=5)
        if show_mutation_label and mut.get("label"):
            ax.text(mut["pos"], 1.05, mut["label"], rotation=45, fontsize=8, ha='left', va='bottom')
    ax.legend(handles=[patches.Patch(color=c, label=t) for t,c in type_colors.items()], bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlim(gene.start-100, gene.end+100)
    ax.set_ylim(-0.2, 1.2)
    ax.set_yticks([])
    ax.set_xlabel(f"{gene.chrom}:{gene.start}-{gene.end}")
    ax.set_title(f"Gene model: {gene.attributes.get('gene_name',[''])[0]}")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    print(f"Gene model figure saved to {output}")

# -----------------------------
# 4. CLI entry
# -----------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Draw gene model with mutations/SVs (CSV or VCF)")
    parser.add_argument("-t","--gtf", required=True, help="GTF/GFF3 gene annotation file")
    parser.add_argument("-g","--gene", required=True, help="Gene symbol or ID to plot")
    parser.add_argument("-m","--mutations", required=False, help="CSV or VCF of mutations")
    parser.add_argument("-o","--output", default="gene_model.png", help="Output figure path")
    parser.add_argument("-f","--figsize", default="10,2", help="Figure size as width,height")
    parser.add_argument("--show_labels", action="store_true", help="Show mutation labels")
    args = parser.parse_args()

    fig_width, fig_height = map(float, args.figsize.split(","))
    db = create_db(args.gtf)
    transcripts, gene = get_gene_structure_by_transcript(db, args.gene)
    mutations = load_mutations(args.mutations)
    outpng = f"{args.output}/{args.gene}_gene_model.png"
    plot_gene_model_all_transcripts(
        gene,
        transcripts,
        mutations=mutations,
        figsize=(fig_width, fig_height),
        output=outpng,
        show_mutation_label=args.show_labels
    )

def run():
    vcfs = {
        "PlaB_P6": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB_only.vcf",
        "PlaB_P20": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB_only.vcf"
    }
    outdirs = {
        "PlaB_P6": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/gene",
        "PlaB_P20": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/gene"
    }
<<<<<<< HEAD
    # genes = ["Fhit","Alk","Lrp1b"]
    genes = ["Ciita","Fgfr2"]
=======
    genes = ["Fhit","Alk","Lrp1b"]
>>>>>>> 0ea0c995979175199379f89ae752af4c876178c6
    db = create_db("/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf")
    for sample, vcf in vcfs.items():
        outdir = outdirs[sample]
        os.makedirs(outdir, exist_ok=True)
        for gene_name in genes:
            outpng = f"{outdir}/{gene_name}_gene_model.png"
            transcripts, gene = get_gene_structure_by_transcript(db, gene_name)
            mutations = load_mutations(vcf)
            outpng = f"{outdir}/{gene_name}_gene_model.png"
            plot_gene_model_all_transcripts(
                gene,
                transcripts,
                mutations=mutations,
                output=outpng,
                show_mutation_label=False
            )

if __name__ == "__main__":
    # main()
    run()