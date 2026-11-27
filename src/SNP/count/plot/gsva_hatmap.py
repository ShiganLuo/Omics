from os import path
from random import sample
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

def plot_heatmap_with_group_colors(
        heatmap_data:pd.DataFrame,
        outplot:str,
        row_rename_map: dict,
        column_rename_map: dict,
        show_group_labels:bool = True,
        title:str = "Heatmap of Positions",
        xlabel:str = "Positions",
        colorbar_width:float = 0.02,
        x_annotation_height:float = 0.2
):
    heatmap_data.rename(index=row_rename_map, inplace=True)
    heatmap_data.rename(columns=column_rename_map, inplace=True)
    heatmap_data = heatmap_data.sort_index()
    new_col_order = [new_name for new_name in column_rename_map.values()]
    print(heatmap_data.head())
    heatmap_data = heatmap_data[new_col_order]
    group_colors = {}
    for group in heatmap_data.index:
        if group not in group_colors:
            group_colors[group] = sample(
                [
                    '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                    '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
                    '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                    '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
                ],
                1
            )[0]
    
    n_rows = heatmap_data.shape[0]
    fig = plt.figure(figsize=(7, 11))
    gs = gridspec.GridSpec(2, 2, width_ratios=[colorbar_width, 1],height_ratios = [1,x_annotation_height],wspace=0,hspace=0.3)
    ax_bar = fig.add_subplot(gs[0, 0])
    group_positions = {}
    for idx, group in enumerate(heatmap_data.index):
        ax_bar.add_patch(
            plt.Rectangle((0.6, idx), 0.3, 1, facecolor=group_colors[group])
        )
        if group not in group_positions:
            group_positions[group] = []
        group_positions[group].append(idx)
    if show_group_labels:
        for group, positions in group_positions.items():
            mid_idx = positions[len(positions)//2]
            ax_bar.text(
                0, mid_idx + 0.5, group,
                ha='right', va='center',
                fontsize=8
            )
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(0, n_rows)
    ax_bar.invert_yaxis()
    ax_bar.axis('off')

    ax_heat = fig.add_subplot(gs[0, 1])
    sns.heatmap(
        heatmap_data,
        ax=ax_heat,
        cmap='coolwarm',
        linewidths=0.01,
    )
    for label in ax_heat.get_xticklabels():
        label.set_rotation(60)
        label.set_fontsize(8)
        label.set_ha('center')
    ax_heat.set_yticks([])
    ax_heat.set_ylabel('')
    ax_heat.set_xlabel(xlabel)
    ax_heat.set_title(title)
    cbar = ax_heat.collections[0].colorbar
    cbar.ax.text(
        0.5, 1.02,
        'ES',
        ha='center', va='bottom',
        rotation=0, fontsize=8,
        transform=cbar.ax.transAxes
    )

    ax_legend = fig.add_subplot(gs[1, :])
    ax_legend.axis('off')
    short_names = list(heatmap_data.columns)
    pathway_rename_map = {v:k for k,v in column_rename_map.items()}
    legend_texts = [
        f"{abbr}: {pathway_rename_map[abbr]}"
        for abbr in short_names
        if abbr in pathway_rename_map
    ]

    max_chars_per_row = 60
    rows = []
    current_row = []
    current_len = 0

    for text in legend_texts:
        text_len = len(text)
        if current_len + text_len > max_chars_per_row and current_row:
            rows.append(current_row)
            current_row = [text]
            current_len = text_len
        else:
            current_row.append(text)
            current_len += text_len + 2
    if current_row:
        rows.append(current_row)
    n_rows = len(rows)
    for row_idx, row in enumerate(rows):
        x_pos = 0
        y_pos = n_rows - row_idx
        for text in row:
            ax_legend.text(
                x_pos, y_pos,
                text,
                ha='left', va='center',
                rotation=0,
                fontsize=8
            )
            x_pos += len(text) + 2

    ax_legend.set_xlim(0, max_chars_per_row)
    ax_legend.set_ylim(0, len(rows))
    plt.savefig(outplot, dpi=300)
    plt.close()



def imshow_heatmap(
        df:pd.DataFrame
):
    fig, ax = plt.subplots()

if __name__ == "__main__":
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/PRJNA663159/count/gsva/gseapy.gene_set.gsva.report.csv")
    df_group = pd.read_csv("/disk5/luosg/Totipotent20251031/data/PRJNA663159/download.csv")
    df_group = df_group.loc[df_group["libraryStrategy"] == "RNA-Seq",:]
    df_group["Sample_name"] = df_group["Sample_name"].str.replace(r"RNA-seq_", "", regex=True)
    sample_rename_map = dict(zip(df_group['Sample_name'],df_group['group']))
    pathway_rename_map = {
        "DNA Double Strand Break Response": "DSBR",
        "HDR through Homologous Recombination (HRR)": "HRR",
        "Nonhomologous End-Joining (NHEJ)": "NHEJ",
        "Mismatch Repair": "MMR",
        "Nucleotide Excision Repair": "NER",
        "Base Excision Repair": "BER",
        "Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template":"TLS",
        "Fanconi Anemia Pathway": "FA",
        "Activation of ATR in response to replication stress": "RS",
        "Telomere Maintenance": "TM",
        "Mitotic Spindle Checkpoint": "MSC",
        "Cell Cycle Checkpoints": "CCC"
    }
    df_wide = df.pivot(index='Name', columns='Term', values='ES')
    print(sample_rename_map)
    plot_heatmap_with_group_colors(df_wide,
                 outplot = "/disk5/luosg/Totipotent20251031/PRJNA663159/count/gsva/heatmap.png",
                 row_rename_map = sample_rename_map,
                 column_rename_map = pathway_rename_map,
                 title="",
                 xlabel="Pathways")
