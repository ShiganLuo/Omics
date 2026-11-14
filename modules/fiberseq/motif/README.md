# fiberseq/motif — FIRE peak 内基序扫描与基序中心谱图子模块

C 形态子模块(两规则:`fiberseq_motif_scan` → `fiberseq_motif_profile`,
加 result 锚点):在 FIRE peak 序列里扫描 MEME 基序,再围绕基序位点汇总
m6A/可及性/核小体谱并计算 footprint 分数。共享父目录环境
`../fiberseq.yaml`,引擎脚本 `modules/fiberseq/bin/motif_scan.py` 与
`motif_profile.py`(经 `ROOT_DIR` 引用)。

## 1. 模块输入输出与参数

### 1.1 输入

元素文件约定见父目录 `modules/fiberseq/README.md`(indir = `common/4_elements`):

| 输入 | 规则 | 说明 |
|---|---|---|
| `{indir}/{sample_id}/{sample_id}.fire_peaks.bed` | scan | 扫描区域 BED |
| `genome.fasta` | scan | 参考基因组 FASTA(需 .fai 索引;`_require` 校验) |
| `genome.motif_file` | scan | MEME 格式基序文件(如 JASPAR 导出;`_require` 校验) |
| `{indir}/{sample_id}/{sample_id}.m6a.bed.gz` | profile | 逐纤维 m6A BED12 |
| `{indir}/{sample_id}/{sample_id}.msp.bed.gz` | profile | 逐纤维 MSP BED12 |
| `{indir}/{sample_id}/{sample_id}.nuc.bed.gz` | profile | 逐纤维核小体 BED12 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.motif_sites.bed               # 基序命中位点 BED
{outdir}/{sample_id}/{sample_id}.motif_profile.tsv              # 基序中心平均谱
{outdir}/{sample_id}/{sample_id}.motif_footprint_scores.tsv     # 逐位点 footprint 分数
{outdir}/{sample_id}/{sample_id}.motif_profile.pdf              # 谱图 PDF
{outdir}/{sample_id}/fiberseq_motif_scan_{sample_id}.sh         # 回执(profile 同目录)
logdir/{sample_id}/fiberseq_motif_scan.log                      # 运行日志(profile.log 同目录)
```

`motif_sites.bed` 列:`chrom start end name score strand`,name 形如
`{motif}@{region_name}`;`motif_profile.tsv` 列:`position fiber_coverage
m6a_frac accessibility_frac nuc_occupancy_frac`(position ∈ [-flank, flank));
`motif_footprint_scores.tsv` 列:`chrom start end motif footprint_score
mean_m6a_frac`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.fiberseq_motif.pwm_threshold` | `-t` | 0.8 | 命中阈值 = 该比例 × PWM 满分 |
| `Params.fiberseq_motif.flank` | `-f` | 250 | 基序中心两侧窗口半宽(bp) |
| `Params.fiberseq_motif.min_fibers` | `-i` | 10 | 位点纳入的最少纤维覆盖 |
| `genome.fasta` | `-f`(scan) | — | 参考 FASTA |
| `genome.motif_file` | `-m`(scan) | — | MEME 基序文件 |
| `samples` | — | [] | 样本列表(result 锚点 expand) |
| `indir` / `outdir` / `logdir` | — | `common/4_elements` / `fiberseq_analysis/motif` / `log` | 目录基路径 |
| `ROOT_DIR` | — | `.` | 仓库根(拼 bin/ 脚本路径) |

线程/内存经 `rule_res("fiberseq_motif")` 配置:scan threads 默认 4
(chunk 并行 `--threads`),mem 4000;profile threads 1,mem 8000。

直接调用示例(引擎脚本):

```bash
python modules/fiberseq/bin/motif_scan.py -m motifs.meme -r sample.fire_peaks.bed \
    -f genome.fa -o out/ --sample s1 -t 0.8 --threads 4
python modules/fiberseq/bin/motif_profile.py -m out/sample.motif_sites.bed \
    -a sample.m6a.bed.gz -s sample.msp.bed.gz -u sample.nuc.bed.gz \
    -o out/ --sample s1 -f 250 -i 10
```

## 2. 算法原理与生物学意义

**motif_scan**:解析 MEME 基序为位置概率矩阵,转 log-odds
PWM(`log2(p/0.25)`,伪计数 0.001),对每条区域序列双链滑窗打分
(N/未知碱基记 −2.0);命中判定 score ≥ `pwm_threshold` × 满分
(满分 = 各位置最大分之和),命中即写 BED 位点。

**motif_profile**:围绕每个位点中心 ± `flank` 窗口,按逐纤维 BED12 汇总
四类逐碱基信号——纤维覆盖(nuc+msp 按纤维并集)、m6A 信号、MSP 可及性、
核小体占用;`min_fibers` 过滤低覆盖位点后输出平均谱与逐位点
footprint 分数:

      footprint = max(左翼 m6A frac) + max(右翼 m6A frac) − Σ(基序中心 m6A frac)

翼侧最大值减中心和:基序被蛋白占据时中心 m6A/可及性信号受"足迹"压制,
分数越高表示该位点越可能被结合。

生物学意义:在 FIRE 调控元件序列里找 TF 基序并看基序中心的单分子信号
剖面,可区分"有基序且被占据"与"有基序但空置"的元件,把序列潜力与真实
染色质状态对应起来。

## 3. 参考文献

1. fiberseq 家族方法来源:见 `modules/fiberseq/README.md` Sources —
   https://doi.org/10.1016/j.crmeth.2024.100911 、
   https://doi.org/10.1101/gr.279095.124 、
   https://doi.org/10.1038/s41588-025-02324-w
2. FIRE analysis pipeline — https://github.com/fiberseq/FIRE
3. MEME motif format / JASPAR — https://meme-suite.org/meme/doc/meme-format.html ;
   https://jaspar.elixir.no/
