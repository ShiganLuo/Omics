# scanpy — 单细胞 RNA-seq 分析模块(QC → 合并 → 聚类注释 → 高级分析 → pseudobulk)

单细胞(scanpy)全流程:逐样本 QC/去双胞,按组织(tissue)合并,自动
聚类 + 注释(支持 LLM 辅助与手工覆写),可选轨迹/RNA velocity/细胞通讯/
CNV 高级分析,最后按细胞类型做 pseudobulk 差异表达(PyDESeq2)并出 HTML
报告。引擎为模块私有脚本 `bin/scRNAseq.py`(qc/merge/auto/reannotate/
advanced 子命令)与 `bin/pseudobulk.py`;注释规范见模块内 `ANNOTATION.md`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 布局参数 | 形态 A 拼接:`{indir}/{sample_id}/{sample_id}_{counter}.h5ad`(scanpy 对象;`counter` 通配符标识计数通道,取值由上游文件名决定) |
| `Params.scanpy.tissue_samples` | 组织 → 样本列表字典(node.py 注入;.smk 读 `Params.scanpy.tissue_samples`,json 示例为顶层 `tissue_samples`,以 .smk 为准) |
| `te_bed` / `gene_tsv` | TE BED 与基因 TSV(合并时做基因类型注释,均给才传) |
| `species` | 物种(组织特异注释用) |
| `Params.scanpy.<counter>.auto.marker_supplement.<tissue>` | 组织特异 marker 文件(可选) |
| `Params.scanpy.<counter>.auto.cluster_annotations.<tissue>` | 手工注释覆写 JSON(启用 scanpy_reannotate) |
| `Params.scanpy.<counter>.pseudobulk.compare` | 比较定义列表 `tissue:ref_sample:treat_sample[=label]`,按组织过滤后传入 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_{counter}_qc.h5ad            # QC 后对象
{outdir}/{sample_id}/{sample_id}_{counter}_qc_metrics.tsv     # QC 指标
{outdir}/{sample_id}/plots/{counter}/                          # QC 图(directory)
{outdir_combine}/{tissue}/{tissue}_{counter}_merged.h5ad      # 组织合并
{outdir_combine}/{tissue}/plots/{counter}/merge/               # 合并图
{outdir_combine}/{tissue}/{tissue}_{counter}_auto.h5ad        # 聚类+注释
{outdir_combine}/{tissue}/{tissue}_{counter}_auto.h5ad.reannotated  # 覆写标记(DAG 排序)
{outdir_combine}/{tissue}/{tissue}_{counter}_advanced.h5ad    # 高级分析
{outdir_combine}/{tissue}/plots/{counter}/advanced/            # 高级分析图
{outdir_pseudobulk}/{tissue}/{counter}/pseudobulk_report.html  # pseudobulk 报告
```

`outdir_combine` 默认 = `outdir`;`outdir_pseudobulk` 默认 = `outdir_combine`
(scRNAseq 约定 `<outdir>/diff/pseudobulk`)。scanpy_reannotate 就地改写
auto.h5ad(备份 `<auto>.bak`,重生成注释报告),只对配置了覆写 JSON 的
组织进 DAG。日志:`{logdir}/{sample_id}/scanpy_qc_{counter}.log`、
`{logdir_combine}/scanpy/{tissue}/scanpy_*_{counter}.log`。

### 1.3 参数

参数按 `Params.scanpy.<counter>.<stage>.<key>` 取值(json 示例为平铺示意):

| 参数 | CLI(scRNAseq.py) | 默认 | 说明 |
|---|---|---|---|
| `qc.min_genes` / `max_genes` | `-n` / `-N` | `200` / `6000` | 基因数上下限 |
| `qc.max_pct_mt` | `-t` | `20` | 线粒体比例上限(%) |
| `qc.n_top_genes` | — | `3000` | 高变基因数 |
| `qc.use_mad` | `-M` | `false` | 用 MAD 自适应阈值 |
| `qc.scrublet` / `doublet_rate` | `-s` / `-d` | `true` / `0.06` | Scrublet 去双胞 |
| `cluster.resolution` | `-r`(auto) | `0.8` | Leiden 分辨率 |
| `auto.max_iterations` | `-I` | `3` | 自动注释迭代上限 |
| `auto.min_genes` / `min_counts` / `max_pct_mt` | `-n` / `-U` / `-t` | `800` / `3000` / `20` | auto 步过滤 |
| `auto.n_pcs` / `n_neighbors` / `n_top_genes` | `-c` / `-k` / `-g` | `50` / `50` / `3000` | 降维/邻居/HVG |
| `auto.batch_method` / `batch_key` | `-B` / `-K` | `harmony` / `sample_id` | 批次校正 |
| `auto.auto_n_pcs` | `-A` | `true` | 自动选主成分数 |
| `auto.skip_te` | `-x` | `true` | 跳过 TE 层分析 |
| `auto.regress_cell_cycle` | `-C`(反向) | `true` | 回归细胞周期(false 时传 `-C`) |
| `auto.llm_method` / `llm_model` / `llm_api_key` / `llm_base_url` | `-L` / `-l` / `-Q` / `-u` | `None` | LLM 辅助注释(见 ANNOTATION.md) |
| `auto.marker_supplement` / `marker_supplement_llm` | `-X` / `-M`(反向) | `{}` / `true` | marker 补充 |
| `auto.exclude_genes` / `exclude_from_deg` | `-E` / `-D` | `""` / `false` | 排除基因 |
| `auto.save_iteration_h5ad` | `-H` | `false` | 保存迭代中间对象 |
| `auto.cluster_annotations.<tissue>` | `-Y`(reannotate) | `[]` | 手工注释覆写 JSON |
| `advanced.trajectory` / `velocity` / `communication` / `cnv` | `-T` / `-V` / `-L` / `-C` | `false`×4(json;.smk 回退 trajectory=true) | 高级分析开关 |
| `advanced.gtf` / `cnv_reference` | `-G` / `-R` | `""` | CNV 注释/参照 |
| `cluster.n_pcs` / `n_neighbors`(advanced) | `-c` / `-k` | `50` / `15` | 高级分析降维参数 |
| `pseudobulk.enabled` | — | `false` | pseudobulk 开关 |
| `pseudobulk.compare` | `--compare`(可重复) | `[]` | `tissue:ref:treat[=label]` |
| `pseudobulk.min_cells` / `min_counts` | `--min-cells` / `--min-counts` | `30` / `1000` | 细胞/计数下限 |
| `pseudobulk.n_pseudo_reps` / `min_cells_per_rep` | `--n-pseudo-reps` / `--min-cells-per-rep` | `3` / `10` | 伪重复设计 |
| `pseudobulk.cell_types` / `top_cell_types` / `highlight_genes` | `--cell-types` / `--top-cell-types` / `--highlight-genes` | `[]` | 目标细胞类型/基因 |
| `pseudobulk.deg_mode` | `--deg-mode`(nargs+) | `["all"]` | DEG 模式 |
| `formats`(qc/auto/advanced/pseudobulk 段内或 counter 级共享) | `--format`(可重复) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集);出图文件 `<base>.<fmt>` 每种格式一个(图目录为 directory(),形态不变) |
| `Procedure.python` | — | `python` | 解释器 |

工程细节:所有 run 块导出 `OMP/OPENBLAS/MKL/NUMBA/VECLIB` 线程上限 = 规则
线程预算;`--cell-types`/`--deg-mode` 是 nargs='+' 单 flag 多值(重复传会
互相覆盖)。环境:`conda: "scanpy.yaml"`(模块目录内)。

直接调用示例(等价引擎命令):

```bash
python modules/scanpy/bin/scRNAseq.py qc -i s1_gene.h5ad -o s1_gene_qc.h5ad \
  -m s1_gene_qc_metrics.tsv -p out/s1/plots/gene -n 200 -N 6000 -t 20 -s -d 0.06
python modules/scanpy/bin/scRNAseq.py auto -i brain_gene_merged.h5ad \
  -o brain_gene_auto.h5ad -T brain -r 0.8 -I 3 -n 800 -U 3000 -t 20 \
  -c 50 -k 50 -g 3000 -B harmony -K sample_id -A -x
python modules/scanpy/bin/pseudobulk.py \
  --h5ad brain=brain_gene_auto.h5ad --counts-h5ad brain=brain_gene_merged.h5ad \
  --out-dir out/diff/pseudobulk/brain/gene \
  --min-cells 30 --min-counts 1000 --n-pseudo-reps 3 --min-cells-per-rep 10 \
  --compare brain:young:old=aged --deg-mode all
```

## 2. 算法原理与生物学意义

- **QC**:过滤基因数/线粒体比例异常细胞,Scrublet(模拟双胞相似度)去
  双胞,是单细胞分析可信度的第一道关。
- **auto**:标准化 → 高变基因 → PCA(可自动选维)→ 邻居图 → Leiden 聚类
  → harmony 批次校正 → 标记基因 + 可选 LLM 辅助细胞类型注释,多轮迭代
  (ANNOTATION.md 为注释规范);`skip_te` 控制是否同时分析 TE 表达层。
- **reannotate**:人工用覆写 JSON 修正聚类标签(就地生效,不重聚类),
  advanced/pseudobulk 通过 `.reannotated` 标记保证读到修正后标签。
- **advanced**:轨迹推断(分化轨迹)、RNA velocity(未来状态)、细胞
  通讯(配体-受体)、CNV(肿瘤拷贝数推断)。
- **pseudobulk**:按细胞类型聚合为伪重复样本跑 PyDESeq2——比逐细胞检验
  更严格控制假阳性(伪重复消除细胞间非独立性);两文件设计:注释取
  auto.h5ad、整数 counts 取 merged.h5ad(auto 的 counts 层是 HVG 子集,
  不能作 pseudobulk 输入)。

结果怎么读:`*_qc_metrics.tsv` 看过滤比例与双胞率;`*_auto.h5ad` 内
leiden 分群与 cell_type 标注(报告目录含 marker 图与注释审计);
`pseudobulk_report.html` 含每细胞类型 DEG(火山/热图/GO)。

## 3. 参考文献

1. scanpy — https://scanpy.readthedocs.io/ ,
   https://github.com/theislab/scanpy(Wolf FA et al., *SCANPY:
   large-scale single-cell gene expression data analysis*. Genome Biology
   19:15, 2018. DOI: 10.1186/s13059-018-1603-1)
2. Scrublet — https://github.com/AllonKleinLab/scrublet(Wolock SL et al.,
   Cell Systems 8:281–291, 2019. DOI: 10.1016/j.cels.2018.11.005)
3. PyDESeq2 — https://github.com/owkin/PyDESeq2
4. Harmony 批次校正 — https://github.com/immunogenomics/harmony
