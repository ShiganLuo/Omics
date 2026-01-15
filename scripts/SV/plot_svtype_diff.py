import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency


def plot_svtype_comparison(
    df,
    out_png,
    group_col="group",
    svtype_col="svtype",
    count_col="count",
    group_order=("Control", "Experiment"),
    svtype_order=("BND", "DEL", "DUP", "INS", "INV"),
    legend_map=None,
    figsize=(9, 5),
    ylabel="SV count",
    dpi=300,
):

    if legend_map is None:
        legend_map = {g: g for g in group_order}

    pivot = (
        df.pivot(index=svtype_col, columns=group_col, values=count_col)
        .reindex(svtype_order)
        .fillna(0)
    )

    def p_to_star(p):
        if p < 1e-4:
            return "****"
        elif p < 1e-3:
            return "***"
        elif p < 1e-2:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"

    stars = {}
    for sv in pivot.index:
        g1, g2 = group_order
        table = [
            [pivot.loc[sv, g1], pivot[g1].sum() - pivot.loc[sv, g1]],
            [pivot.loc[sv, g2], pivot[g2].sum() - pivot.loc[sv, g2]],
        ]
        _, p, _, _ = chi2_contingency(table)
        stars[sv] = p_to_star(p)

    sig_sv = [sv for sv, s in stars.items() if s != "ns"]
    low_max = (
        pivot.loc[sig_sv].values.max() * 1.15
        if sig_sv else
        np.median(pivot.values) * 1.5
    )

    global_max = pivot.values.max()
    high_min = low_max * 1.1
    high_max = global_max * 1.15   # ⬅️ 预留空间给显著帽子

    x = np.arange(len(pivot.index))
    width = 0.36

    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1, sharex=True,
        figsize=figsize,
        gridspec_kw={"height_ratios": [1, 3]},
    )

    colors = {
        group_order[0]: "#4C72B0",
        group_order[1]: "#DD8452",
    }

    for ax in (ax_top, ax_bottom):
        ax.bar(x - width / 2, pivot[group_order[0]], width,
               color=colors[group_order[0]], label=legend_map[group_order[0]])
        ax.bar(x + width / 2, pivot[group_order[1]], width,
               color=colors[group_order[1]], label=legend_map[group_order[1]])

    ax_bottom.set_ylim(0, low_max)
    ax_top.set_ylim(high_min, high_max)

    # ---------- fixed broken axis (LEFT ONLY, SAFE) ----------
    d = 0.008

    # top panel: bottom edge
    ax_top.plot(
        (-d, +d),
        (-d, +d),
        transform=ax_top.transAxes,
        color="black",
        clip_on=False,
    )

    # bottom panel: top edge
    ax_bottom.plot(
        (-d, +d),
        (1 - d, 1 + d),
        transform=ax_bottom.transAxes,
        color="black",
        clip_on=False,
    )


    # ---------- significance (VISUALLY CONSISTENT, SAFE) ----------
    LEG_PT = 8     # 显著腿高度（物理单位）
    TEXT_PT = 3

    fig.canvas.draw()  # 保证 transform 可用

    for i, sv in enumerate(pivot.index):
        y1 = pivot.loc[sv, group_order[0]]
        y2 = pivot.loc[sv, group_order[1]]
        y_base = max(y1, y2)

        ax = ax_bottom if y_base <= low_max else ax_top

        # data -> display
        trans = ax.transData
        inv = ax.transData.inverted()

        _, y_disp = trans.transform((0, y_base))
        _, y_hat = inv.transform((0, y_disp + LEG_PT))
        _, y_text = inv.transform((0, y_disp + LEG_PT + TEXT_PT))

        x1 = x[i] - width / 2
        x2 = x[i] + width / 2

        ax.plot([x1, x1], [y1, y_hat], lw=1.2, c="black")
        ax.plot([x2, x2], [y2, y_hat], lw=1.2, c="black")
        ax.plot([x1, x2], [y_hat, y_hat], lw=1.2, c="black")

        ax.text(
            x[i], y_text, stars[sv],
            ha="center", va="bottom",
            fontsize=12,
            fontweight="bold",
        )

    ax_bottom.set_xticks(x)
    ax_bottom.set_xticklabels(pivot.index)
    ax_bottom.set_ylabel(ylabel)

    ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_top.spines["bottom"].set_visible(False)

    for ax in (ax_top, ax_bottom):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax_top.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.close()


if __name__ == "__main__":
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/PacBio/DEG/table/sv_type_summary.tsv",sep="\t")
    plot_svtype_comparison(
        df,
        "/disk5/luosg/Totipotent20251031/PacBio/DEG/plot/sv_type.png",
        legend_map={
        "Control": "DMSO",
        "Experiment": "PlaB"
    })
