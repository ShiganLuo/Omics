# spatial_scanpy — 空间转录组 Scanpy 分析模块

对 Visium/自定义空间转录组数据做 QC → 聚类 → 空间自相关三段分析
(`spatial_scanpy_qc` → `spatial_scanpy_cluster` → `spatial_scanpy_advanced`),
引擎为模块自带 `bin/spatial_transcriptomics.py`(scanpy + squidpy)。

## 1. 模块输入输出与参数

### 1.1 输入

二选一(config 顶层键;两者皆空时模块加载即抛 `ValueError`):

| 输入 | 说明 |
|---|---|
| `input_h5ad` | 已有 AnnData `.h5ad`,必须含 `adata.obsm["spatial"]` 坐标(脚本校验) |
| `visium_h5` + `spatial_dir` | 10x Visium `.h5` 计数文件;`spatial_dir` 非空时用 `sc.read_visium` 连图像一起读(spatial_dir 为空则 `sc.read_10x_h5`) |

### 1.2 输出

```
{outdir}/spatial/qc/filtered.h5ad                        # QC 后 AnnData(留 counts 层)
{outdir}/spatial/qc/qc_metrics.tsv                       # QC 指标表(n_genes/total_counts/pct_mt)
{outdir}/spatial/cluster/clustered.h5ad                  # 聚类后 AnnData(spatial_leiden 注释)
{outdir}/spatial/cluster/spatial_markers.tsv             # Leiden 分群 marker 基因(Wilcoxon)
{outdir}/spatial/advanced/spatial_autocorrelation.h5ad   # Moran's I 结果(进 obsm/var)
{outdir}/spatial/qc/spatial_qc.sh 等                     # 各规则回执
logdir/spatial/spatial_scanpy_qc.log 等                  # 各规则日志
```

`spatial_scanpy_advanced` 仅当 `Params.spatial.advanced.spatial_autocorrelation`
非 false 时注册(.smk 顶层条件)。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.python` | 可执行文件 | `python` | 解释器 |
| `input_h5ad` | `--input-h5ad` | "" | 输入 h5ad |
| `visium_h5` / `spatial_dir` | `--visium-h5` / `--spatial-dir` | "" | Visium 输入 |
| `Params.spatial.min_genes` | `--min-genes` | 100 | spot 最少基因数(qc) |
| `Params.spatial.max_pct_mt` | `--max-pct-mt` | 25 | 线粒体比例上限 %(qc) |
| `Params.spatial.n_top_genes` | `--n-top-genes` | 3000 | 高变基因数(qc) |
| `Params.spatial.n_pcs` | `--n-pcs` | 50 | PCA 维数(cluster) |
| `Params.spatial.n_neighbors` | `--n-neighbors` | 15 | 近邻数(cluster) |
| `Params.spatial.resolution` | `--resolution` | 0.8 | Leiden 分辨率(cluster) |
| `Params.spatial.advanced.spatial_autocorrelation` | — | true | 是否注册 advanced 规则 |
| `Params.spatial.advanced.deconvolution` | — | false | json 预留键,.smk 未读取 |
| `ROOT_DIR` / `outdir` / `logdir` | — | `.` / `output` / `logs` | 基础键 |

线程/内存经 `rule_res(...)` 配置:qc/cluster threads 4、mem 8000
(advanced 固定 4 线程)。

直接调用示例(模块 bin/ 引擎脚本):

```bash
python modules/spatial_scanpy/bin/spatial_transcriptomics.py --mode qc \
  --visium-h5 raw_feature_bc_matrix.h5 --spatial-dir spatial/ \
  --output out/spatial/qc/filtered.h5ad --metrics out/spatial/qc/qc_metrics.tsv \
  --min-genes 100 --max-pct-mt 25 --n-top-genes 3000
python modules/spatial_scanpy/bin/spatial_transcriptomics.py --mode cluster \
  --input-h5ad out/spatial/qc/filtered.h5ad --output out/spatial/cluster/clustered.h5ad \
  --genes out/spatial/cluster/spatial_markers.tsv --n-pcs 50 --n-neighbors 15 --resolution 0.8
```

## 2. 算法原理与生物学意义

- **qc**:标 MT- 基因为线粒体,`calculate_qc_metrics` 后过滤
  (`n_genes_by_counts ≥ min_genes`、`pct_counts_mt ≤ max_pct_mt`),
  保存 counts 层,`normalize_total(1e4)` + `log1p`,挑高变基因;
  qc_metrics.tsv 供人工核查过滤效果。
- **cluster**:取高变基因 → `scale(max_value=10)` → PCA → 近邻图 →
  UMAP → Leiden 聚类(`spatial_leiden`),Wilcoxon 秩和检验输出各群
  marker(spatial_markers.tsv)。
- **advanced**:squidpy 建空间近邻图,对前 2000 个基因算 Moran's I
  空间自相关——识别表达随组织空间位置系统变化的基因(空间高变基因)。

生物学意义:空间转录组把表达量锚定在组织切片坐标上;QC 去除低质量
spot 与死细胞,Leiden 分群对应组织学区域,Moran's I 把"表达随位置变化"
的基因挑出来,用于空间域标志物与微环境分析。

## 3. 参考文献

1. Scanpy — https://scanpy.readthedocs.io/ ;Wolf FA et al. *Scanpy:
   large-scale single-cell gene expression data analysis*. Genome Biology
   19:15, 2018. DOI: 10.1186/s13059-017-1382-0
2. Squidpy — https://github.com/scverse/squidpy ;Palla G et al. *Squidpy:
   a scalable framework for spatial omics analysis*. Nature Methods
   19(2):171–178, 2022. DOI: 10.1038/s41592-021-01358-2
