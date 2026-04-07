#!/usr/bin/env python3
# TE_SNV_full_pipeline_parallel_stream.py
# Parallel + streaming, low-memory version of TE-SNV pipeline

import os
import argparse
import gzip
import multiprocessing as mp
from functools import partial
import tempfile
import shutil

import pandas as pd
import numpy as np
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ---------------------------
# Styling
# ---------------------------
def set_nature_style():
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8,
    })

# ---------------------------
# Natural chrom sort key
# ---------------------------
def chrom_key(ch):
    c = str(ch)
    if c.startswith("chr"):
        name = c[3:]
    else:
        name = c
    if name.isdigit():
        return (0, int(name))
    if name in ("X","x"):
        return (1, 23)
    if name in ("Y","y"):
        return (1, 24)
    if name in ("M","MT","m","Mt"):
        return (1, 25)
    return (2, name)

def sort_df_by_chrom(df, chrom_col="chrom", start_col="start"):
    df2 = df.copy()
    df2["_ck"] = df2[chrom_col].map(chrom_key)
    df2 = df2.sort_values(["_ck", start_col])
    return df2.drop(columns=["_ck"])

# ---------------------------
# Load TE GTF -> BedTool (4 columns: chrom,start,end,attr)
# This returns a BedTool with proper sorting.
# ---------------------------
def load_te_gtf_as_bed_sorted(gtf_path, tmp_dir=None):
    rows = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, _, _, attr = parts[:9]
            try:
                s = int(start) - 1
                e = int(end)
            except:
                continue
            rows.append([chrom, s, e, attr])
    if not rows:
        raise SystemExit("No TE records parsed from GTF.")
    df = pd.DataFrame(rows, columns=["chrom","start","end","attr"])
    df = sort_df_by_chrom(df, "chrom", "start")
    # create a temp file and return BedTool from dataframe
    if tmp_dir:
        tmpfn = os.path.join(tmp_dir, "te_sorted.bed")
        df.to_csv(tmpfn, sep="\t", index=False, header=False)
        bt = pybedtools.BedTool(tmpfn)
    else:
        bt = pybedtools.BedTool.from_dataframe(df)
    # final sort via bedtools, return object
    return bt.sort()

# ---------------------------
# Load a VCF into a sorted BedTool (streamed)
# Returns a BedTool; creates temporary file on disk to keep memory low.
# ---------------------------
def vcf_to_bed_sorted_file(vcf_path, tmp_dir):
    # Write a temp BED-like file with columns chrom start end ref alt
    basename = os.path.basename(vcf_path).replace(".vcf.gz","").replace(".vcf","")
    tmpfn = os.path.join(tmp_dir, f"{basename}.bed")
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh, open(tmpfn, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except:
                continue
            ref = parts[3]
            alt = parts[4]
            start = pos - 1
            end = pos
            out.write(f"{chrom}\t{start}\t{end}\t{ref}\t{alt}\n")
    # sort using pybedtools/bedtools (writes a sorted temp file)
    bt = pybedtools.BedTool(tmpfn).sort()
    return bt

# ---------------------------
# Robust closest wrapper per-sample
# This will run in parallel in worker processes.
# It writes per-sample TSV directly to disk (streaming), returns path.
# ---------------------------
def process_sample_closest(sample, vcf_path, te_bed, outdir, tmp_dir):
    """
    sample: sample name
    vcf_path: path to vcf file
    te_bed: pybedtools.BedTool (can be pickled? pybedtools objects are not picklable)
            therefore we instead require te_bed_path (filepath) OR we build te_bed anew in each worker
    outdir: directory to write sample output
    tmp_dir: directory for temporary files
    """
    # Because BedTool objects are not picklable, this function will accept te_bed as a filepath string
    # In our caller we'll pass te_bed_path rather than BedTool object.
    raise RuntimeError("do not call directly; use worker wrapper below")

# Worker wrapper that accepts te_bed_path
def worker_process_sample(args):
    sample, vcf_path, te_bed_path, outdir, tmp_dir = args
    try:
        # prepare vcf bed sorted
        vcf_bt = vcf_to_bed_sorted_file(vcf_path, tmp_dir)
        # prepare te bed sorted from file (each worker makes its own BedTool from file path; cheap)
        te_bt = pybedtools.BedTool(te_bed_path).sort()
        # run bedtools closest -d
        closest = vcf_bt.closest(te_bt, d=True).to_dataframe(disable_auto_names=True)
        # robust column parsing: last column is distance
        cols = closest.columns.tolist()
        dist_col = cols[-1]
        # Map vcf first 5 columns if available
        mapping = {}
        if len(cols) >= 5:
            mapping[cols[0]] = "v_chrom"
            mapping[cols[1]] = "v_start"
            mapping[cols[2]] = "v_end"
            mapping[cols[3]] = "ref"
            mapping[cols[4]] = "alt"
        mapping[dist_col] = "distance"
        # remaining -> te columns
        te_cols = [c for c in cols if c not in mapping]
        # name te columns
        te_names = ["t_chrom","t_start","t_end","t_attr"]
        for i, c in enumerate(te_cols):
            mapping[c] = te_names[i] if i < len(te_names) else f"t_col_{i}"
        closest = closest.rename(columns=mapping)
        # ensure columns exist
        for c in ["v_chrom","v_start","v_end","ref","alt","t_chrom","t_start","t_end","t_attr","distance"]:
            if c not in closest.columns:
                closest[c] = np.nan
        # enforce numeric
        closest["v_start"] = pd.to_numeric(closest["v_start"], errors="coerce").fillna(-1).astype(int)
        closest["v_end"] = pd.to_numeric(closest["v_end"], errors="coerce").fillna(-1).astype(int)
        closest["t_start"] = pd.to_numeric(closest["t_start"], errors="coerce").fillna(-1).astype(int)
        closest["t_end"] = pd.to_numeric(closest["t_end"], errors="coerce").fillna(-1).astype(int)
        closest["distance"] = pd.to_numeric(closest["distance"], errors="coerce").fillna(-1).astype(int)
        # set inside TE -> distance 0
        mask_inside = (
            (closest["v_chrom"] == closest["t_chrom"]) &
            (closest["v_start"] >= closest["t_start"]) &
            (closest["v_end"] <= closest["t_end"])
        )
        closest.loc[mask_inside, "distance"] = 0
        # add sample name
        closest["sample"] = sample
        # write per-sample file (tsv)
        outfn = os.path.join(outdir, f"{sample}.distance.tsv")
        # keep consistent columns order
        cols_out = ["sample","v_chrom","v_start","v_end","ref","alt","t_chrom","t_start","t_end","t_attr","distance"]
        # ensure columns present
        for c in cols_out:
            if c not in closest.columns:
                closest[c] = ""
        closest[cols_out].to_csv(outfn, sep="\t", index=False)
        # cleanup worker temp if any (vcf_bt uses tmp file created in tmp_dir)
        return (sample, outfn, None)
    except Exception as e:
        return (sample, None, str(e))

# ---------------------------
# Stream-merge per-sample tsvs into one big TSV without loading all into memory
# ---------------------------
def stream_merge_tsvs(sample_files, merged_outpath, chunk_lines=100000):
    """
    sample_files: list of file paths
    merged_outpath: destination file
    chunk_lines: lines to read per chunk (for pandas read_csv chunksize)
    """
    first = True
    # We'll iterate over files and write in append mode
    with open(merged_outpath, "w") as outfh:
        for i, f in enumerate(sample_files):
            # read header to skip after first
            with open(f, "r") as infh:
                header = infh.readline().rstrip("\n")
                if first:
                    outfh.write(header + "\n")
                    first = False
                # stream rest lines in chunks
                # simply copy lines (avoid pandas overhead)
                for line in infh:
                    outfh.write(line)
    return merged_outpath

# ---------------------------
# Plotting helpers (memory-conscious)
# For plotting we will sample a fraction of merged file if it's huge.
# ---------------------------
def load_merged_sample(merged_path, sample_frac=0.2, seed=42, max_rows=None):
    """
    Returns a DataFrame loaded by sampling lines from the TSV.
    For large files, sampling avoids loading entire file.
    If sample_frac == 1.0, will try to load whole file (may be memory heavy).
    """
    if sample_frac >= 1.0:
        # load full (optionally cap by max_rows)
        df = pd.read_csv(merged_path, sep="\t")
        if max_rows and len(df) > max_rows:
            return df.sample(n=max_rows, random_state=seed)
        return df
    # streaming sampling: read line by line and keep a reservoir sample
    import random
    random.seed(seed)
    sample = []
    header = None
    total = 0
    with open(merged_path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            total += 1
            if random.random() < sample_frac:
                sample.append(line.rstrip("\n").split("\t"))
            # optional cap
            if max_rows and len(sample) >= max_rows:
                break
    if not sample:
        # fallback: read small chunk
        df = pd.read_csv(merged_path, sep="\t", nrows=10000)
        return df
    df = pd.DataFrame(sample, columns=header)
    return df

# plotting functions (reuse from previous version, but operate on sample df)
def plot_te_internal_box_from_df(df, outpng):
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce").fillna(-1).astype(int)
    df["in_te"] = (df["distance"] == 0).astype(int)
    sample_counts = df.groupby(["sample","group"], as_index=False)["in_te"].sum()
    plt.figure(figsize=(6,4))
    sns.boxplot(data=sample_counts, x="group", y="in_te", color="white")
    sns.stripplot(data=sample_counts, x="group", y="in_te", color=".25", jitter=True, size=4)
    plt.ylabel("SNVs inside TE (count per sample)")
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def plot_distance_kde_from_df(df, outpng):
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce").fillna(-1).astype(int)
    df = df[df["distance"] >= 0]
    if df.empty:
        return
    cap = int(df["distance"].quantile(0.99))
    df["distance_capped"] = df["distance"].clip(upper=cap)
    plt.figure(figsize=(6,4))
    sns.kdeplot(data=df, x="distance_capped", hue="group", common_norm=False, bw_adjust=0.6)
    plt.xlim(0, cap)
    plt.tight_layout()
    plt.savefig(outpng,dpi=300)
    plt.close()

def plot_linear_log_hist_from_df(df, out_linear, out_log):
    df = df.copy()
    vals = df["distance"].astype(float)
    # 合法距离必须 >= 0
    vals = vals[(vals >= 0) & np.isfinite(vals)]

    if len(vals) == 0:
        print("No valid distance values for plotting.")
        return
    sns.set_theme(style="white", font_scale=1.4)

    plt.figure(figsize=(6,4))
    plt.hist(vals, bins=100, edgecolor="black")
    plt.xlabel("Distance to nearest TE (bp)")
    plt.ylabel("Count")
    sns.despine()
    plt.tight_layout()
    plt.savefig(out_linear, dpi=300)
    plt.close()

    log_vals = np.log10(vals + 1)
    log_vals = log_vals[np.isfinite(log_vals)]

    plt.figure(figsize=(6,4))
    plt.hist(log_vals, bins=100, edgecolor="black")
    plt.xlabel("log10(Distance + 1)")
    plt.ylabel("Count")
    sns.despine()
    plt.tight_layout()
    plt.savefig(out_log, dpi=300)
    plt.close()

def plot_family_enrichment_from_df(df: pd.DataFrame, outpng: str, top_n: int = 20):
    """
    绘制所有样本中 TE 内 SNV 前 top_n 个 family 的富集条形图。

    参数：
        df (pd.DataFrame): 包含 "distance" 和 "t_attr" 列的数据。
        outpng (str): 输出图片路径。
        top_n (int): 显示前 top_n 个 family。
    """
    # 提取 family_id
    fams = (
        df.loc[df["distance"] == 0, "t_attr"]
        .astype(str)
        .str.extract(r'family_id "([^"]+)"')[0]
        .fillna("NA")
    )

    # 统计 top_n family
    top = fams.value_counts().head(top_n)
    if top.empty:
        return

    res = top.reset_index()
    res.columns = ["family", "count"]

    # 设置图形大小（高度随条数自动调整）
    plt.figure(figsize=(6, max(3, len(res) * 0.4)))

    # 使用渐变色条形图，更美观
    palette = sns.color_palette("viridis", n_colors=len(res))
    sns.barplot(data=res, y="family", x="count", palette=palette)

    # 美化图表
    plt.xlabel("Count", fontsize=12)
    plt.ylabel("Family", fontsize=12)
    # plt.title(f"Top {top_n} Families in TE-overlapping SNVs", fontsize=14)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()

    # 保存并关闭
    plt.savefig(outpng, dpi=300)
    plt.close()

def plot_manhattan_from_df(df, outpng):
    if not {"v_chrom","v_start","distance"}.issubset(df.columns):
        return
    main_chroms = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"] 
    df2 = df[df["v_chrom"].isin(main_chroms)].copy()
    df2["v_start"] = pd.to_numeric(df2["v_start"], errors="coerce").fillna(0).astype(int)
    chroms = list(dict.fromkeys(df2["v_chrom"].tolist()))
    chroms_sorted = sorted(chroms, key=chrom_key)
    chrom_to_offset = {}
    offset = 0
    ticks=[]
    tickpos=[]
    for ch in chroms_sorted:
        sub = df2[df2["v_chrom"]==ch]
        if sub.empty: continue
        minp = sub["v_start"].min()
        maxp = sub["v_start"].max()
        length = maxp - minp + 1
        chrom_to_offset[ch] = offset
        ticks.append(ch)
        tickpos.append(offset + length/2)
        offset += length + 1000
    df2["gpos"] = df2.apply(lambda r: chrom_to_offset.get(r["v_chrom"],0) + int(r["v_start"]), axis=1)
    plt.figure(figsize=(12,3))
    plt.scatter(df2["gpos"], df2["distance"], s=2, alpha=0.6)
    plt.ylim(0, np.percentile(df2["distance"].values, 99))
    plt.xticks(tickpos, ticks, rotation=90, fontsize=7)
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def plot_family_group_heatmap_from_df(df, outpng, top_n=10):
    """
    绘制 TE 内 SNV family × group 热图，显示归一化比例
    y 轴 family 按总比例排序，使颜色从上到下依次变浅
    """
    # 提取 family
    df["family"] = (
        df["t_attr"].astype(str)
        .str.extract(r'family_id "([^"]+)"')[0]
        .fillna("NA")
    )

    # 只从 TE 内 SNV (distance==0) 中取 top_n family
    top = (
        df[df["distance"] == 0]["family"]
        .value_counts()
        .head(top_n)
        .index
        .tolist()
    )
    if not top:
        return

    # 子集仅包含 TE 内 SNV 且在 top family
    sub = df[(df["distance"] == 0) & (df["family"].isin(top))]

    # group × family 的计数矩阵
    pivot = sub.groupby(["family", "group"]).size().unstack(fill_value=0)

    # 每列归一化
    pivot_norm = pivot.div(pivot.sum(axis=0).replace(0, 1), axis=1)

    # 按总比例排序 family（行）  
    row_order = pivot_norm.sum(axis=1).sort_values(ascending=False).index
    pivot_norm = pivot_norm.loc[row_order]

    # 动态 figsize
    figsize = (
        max(6, len(pivot_norm.columns) * 0.8),
        max(4, len(pivot_norm.index) * 0.5)
    )
    plt.figure(figsize=figsize)
    vmax = pivot_norm.max().max()

    # 绘制 heatmap
    sns.heatmap(
        pivot_norm,
        cmap="rocket_r",        # 颜色对比明显
        cbar_kws={"label": "Proportion"},
        vmin=0,
        vmax=vmax,
        yticklabels=True,
        linewidths=0.5,
        linecolor="gray",
        annot=False
    )

    # 保证 y 轴显示所有标签，字体可调
    plt.yticks(rotation=0, fontsize=8)
    plt.xticks(rotation=45, fontsize=8)

    plt.xlabel("Group",fontsize=10)
    plt.ylabel("Family",fontsize=10)

    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()


def plot_radar_samples_from_df(df, outdir, top_n=8):
    """
    每个样本在 TE 内的 SNV 分布按 TE family 的贡献度, 展示top_n个family
    """
    df["family"] = df["t_attr"].astype(str).str.extract(r'family_id "([^"]+)"')[0].fillna("NA")
    top = df[df["distance"]==0]["family"].value_counts().head(top_n).index.tolist()
    if not top: return
    import math
    for sample, sdf in df.groupby("sample"):
        total_inside = (sdf["distance"]==0).sum()
        if total_inside==0: continue
        vals = [((sdf["distance"]==0) & (sdf["family"]==f)).sum()/total_inside for f in top]
        angles = np.linspace(0, 2*math.pi, len(top), endpoint=False).tolist()
        vals2 = vals + vals[:1]; angles2 = angles + angles[:1]
        plt.figure(figsize=(5,5))
        ax = plt.subplot(111, polar=True)
        ax.plot(angles2, vals2, 'o-', linewidth=1)
        ax.fill(angles2, vals2, alpha=0.25)
        ax.set_thetagrids(np.degrees(angles), top)
        ax.set_ylim(0, max(0.1, max(vals)*1.2))
        outpng = os.path.join(outdir, f"{sample}_radar.png")
        plt.tight_layout()
        plt.savefig(outpng)
        plt.close()



def plot_family_group_barstack_from_df(
    df: pd.DataFrame,
    outpng: str,
    top_n: int = 10,
    figsize: tuple = (10, 5),
    cmap: str = "tab20",
    dpi: int = 300
):
    """
    绘制不同 group 中 TE 内 (distance==0) SNV 的前 top_n 个 family 的堆叠条形图。
    修复了之前的 bug：确保统计与绘图均仅基于 distance == 0 的记录。

    参数：
        df: 包含至少列 ["t_attr","distance","group"] 的 DataFrame
        outpng: 输出图片路径
        top_n: 选取前 top_n 个 family（按 TE 内 SNV 总数排序）
        figsize: 图像大小
        cmap: seaborn/matplotlib 调色板名称
        dpi: 输出图片分辨率
    """
    # 复制输入避免修改原始 df
    dfc = df.copy()

    # 提取 family 字段（保持与原来行为一致）
    dfc["family"] = (
        dfc["t_attr"].astype(str)
        .str.extract(r'family_id "([^"]+)"')[0]
        .fillna("NA")
    )

    # 只在 TE 内 (distance == 0) 统计 top_n family
    te_mask = dfc["distance"] == 0
    fam_counts = dfc.loc[te_mask, "family"].value_counts().head(top_n)
    if fam_counts.empty:
        # 没有 TE 内的 SNV
        return

    top_families = fam_counts.index.tolist()

    # 对绘图数据也只保留 TE 内 (distance==0) 并且属于 top_families 的记录
    sub = dfc.loc[te_mask & dfc["family"].isin(top_families), ["group", "family"]]

    # group x family 计数透视表，确保列顺序与 top_families 保持一致
    tab = sub.groupby(["group", "family"]).size().unstack(fill_value=0)
    # reindex 保证即使某些 family 在某些 group 中全为 0，列顺序仍与 top_families 一致
    tab = tab.reindex(columns=top_families, fill_value=0)

    # 若无 group（例如 group 列缺失或全部 NaN），将使用单一组名
    if tab.shape[0] == 0:
        return

    # 颜色与列严格对应
    colors = sns.color_palette(cmap, n_colors=len(tab.columns))

    ax = tab.plot(
        kind="bar",
        stacked=True,
        figsize=figsize,
        width=0.8,
        color=colors
    )

    # 获取 handles 并与 columns 一一对应，手动设置 legend 保证一致性
    handles, _ = ax.get_legend_handles_labels()
    # pandas/matplotlib 堆叠绘图通常会将第一个列画在底部，legend 顺序与绘制顺序一致，
    # 为了在 legend 中与颜色直观对应，这里使用 handles 与 tab.columns 对齐
    ax.legend(handles, tab.columns, title="Family",
              bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10,frameon = False)

    ax.set_xlabel("Group", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    # ax.set_title(f"Top {len(top_families)} Families in TE-overlapping SNVs by Group", fontsize=14)

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(outpng, dpi=dpi, bbox_inches='tight')
    plt.close()



# ---------------------------
# Main driver
# ---------------------------
def run_pipeline(vcf_dir, group_file, te_gtf, outdir, nproc=4, sample_frac=0.2, stream_chunk=100000):
    os.makedirs(outdir, exist_ok=True)
    tmp_dir = tempfile.mkdtemp(prefix="te_snv_tmp_")
    print("Temporary dir:", tmp_dir)
    # set_nature_style()

    groups = pd.read_csv(group_file, sep="\t", dtype=str)
    if not {"sample","group"}.issubset(groups.columns):
        raise SystemExit("group file must have header sample<TAB>group")

    # prepare sorted TE Bed file on disk to be used by workers
    te_bt = load_te_gtf_as_bed_sorted(te_gtf, tmp_dir=tmp_dir)
    te_bed_path = os.path.join(tmp_dir, "te_sorted.bed")
    # ensure te_bt saved to file (pybedtools writes to temp file location accessible)
    te_bt.saveas(te_bed_path)

    # build worker args list
    worker_args = []
    sample_files = []
    for _, row in groups.iterrows():
        sample = row["sample"]
        group = row["group"]
        # expect sample.vcf.gz (case A)
        p_gz = os.path.join(vcf_dir, f"{sample}.vcf.gz")
        p_plain = os.path.join(vcf_dir, f"{sample}.vcf")
        vpath = p_gz if os.path.exists(p_gz) else (p_plain if os.path.exists(p_plain) else None)
        if vpath is None:
            print(f"[WARN] VCF not found for {sample}, skipping.")
            continue
        out_sample = os.path.join(outdir, f"{sample}.distance.tsv")
        worker_args.append((sample, vpath, te_bed_path, outdir, tmp_dir))
        sample_files.append(out_sample)

    if not worker_args:
        raise SystemExit("No samples found to process.")

    # parallel map
    print(f"Running {len(worker_args)} samples with {nproc} processes...")
    with mp.Pool(processes=nproc) as pool:
        results = pool.map(worker_process_sample, worker_args)

    # report results and collect produced files in same order as groups
    produced = []
    errors = []
    for sample, outfn, err in results:
        if outfn:
            produced.append(outfn)
        else:
            errors.append((sample, err))
    if errors:
        print("Some samples failed:")
        for s,e in errors:
            print(" -", s, e)

    # stream-merge produced files into merged.tsv
    merged_out = os.path.join(outdir, "All_samples.distance.tsv")
    stream_merge_tsvs(produced, merged_out, chunk_lines=stream_chunk)
    print("Merged TSV written to:", merged_out)

    # load sampled subset for plotting / stats (memory friendly)
    print("Sampling merged file for plotting...")
    df_sample = load_merged_sample(merged_out, sample_frac=sample_frac, max_rows=1000000)  # cap to 1e6 rows
    # ensure group column exists
    if "group" not in df_sample.columns:
        # maybe the per-sample files didn't include group; add by parsing filename
        # but earlier worker wrote sample column; we'll attach group via groups map
        sample2group = dict(zip(groups["sample"], groups["group"]))
        if "sample" in df_sample.columns:
            df_sample["group"] = df_sample["sample"].map(sample2group).fillna("unknown")
    # plotting
    plot_te_internal_box_from_df(df_sample, os.path.join(outdir, "te_internal_boxplot.png"))
    plot_distance_kde_from_df(df_sample, os.path.join(outdir, "distance_kde_by_group.png"))
    plot_linear_log_hist_from_df(df_sample, os.path.join(outdir, "distance_linear.png"), os.path.join(outdir,"distance_log10.png"))
    plot_family_enrichment_from_df(df_sample, os.path.join(outdir, "family_enrichment_top20.png"))
    plot_manhattan_from_df(df_sample, os.path.join(outdir, "manhattan_te_distance.png"))
    plot_family_group_heatmap_from_df(df_sample, os.path.join(outdir, "family_group_heatmap.png"))
    plot_radar_samples_from_df(df_sample, outdir, top_n=8)
    plot_family_group_barstack_from_df(df_sample, os.path.join(outdir, "family_group_stack.png"))

    # save sample-level summary (streamed compute to avoid big df)
    print("Computing sample summary (streamed)...")
    # We'll compute per-sample counts by reading produced files sequentially
    sample_stats = []
    for f in produced:
        df_chunk = pd.read_csv(f, sep="\t", usecols=["sample","distance"])
        n_total = len(df_chunk)
        n_in = int((df_chunk["distance"]==0).sum())
        sample_name = os.path.basename(f).split(".")[0]
        sample_stats.append({"sample": sample_name, "n_variants": n_total, "n_in_te": n_in, "pct_in_te": n_in / n_total if n_total>0 else 0})
    pd.DataFrame(sample_stats).to_csv(os.path.join(outdir, "sample_summary_stats.tsv"), sep="\t", index=False)

    # cleanup tmp dir
    try:
        shutil.rmtree(tmp_dir)
    except Exception:
        pass

    print("All done. Outputs in:", outdir)
    if errors:
        print("[WARN] Errors encountered for some samples (see above).")

# ---------------------------
# CLI
# ---------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TE-SNV pipeline: parallel + streaming")
    parser.add_argument("--vcf_dir", required=True)
    parser.add_argument("--group_file", required=True)
    parser.add_argument("--te_gtf", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--nproc", type=int, default=4, help="number of parallel processes")
    parser.add_argument("--sample_frac", type=float, default=0.2, help="fraction of merged file to sample for plotting (0-1)")
    parser.add_argument("--stream_chunk", type=int, default=100000, help="chunk lines for stream merge (unused but kept)")
    args = parser.parse_args()
    run_pipeline(args.vcf_dir, args.group_file, args.te_gtf, args.outdir, args.nproc, args.sample_frac, args.stream_chunk)
