# pseudobulk_go — scRNAseq pseudobulk 结果的 GO/KEGG 富集子模块

`modules/function` 的布局 C 子工具。扫描 scRNAseq pseudobulk 输出树，对每个
`(counter, mode, tissue, contrast, cell_type)` 单元的 DEG 表做 up/down 分开的
GO/KEGG 富集，引擎复用父目录 `modules/function/bin/go-kegg.r`（clusterProfiler +
org.Hs/Mm/Mmu.eg.db，human/mouse/macaque）。te 模式同样受理（TE 名称无法映射
基因库时产空表，属正常结果）。

## 1. 模块输入输出与参数

### 1.1 输入（pseudobulk 输出契约）

```
<root>/<tissue>/<counter>/deg_pseudobulk{,_gene,_te}/
    deg_<tissue>/<contrast>/<cell_type>/<contrast>_<cell_type>_DEG.csv
```

`root` 空则取 `indir`（subworkflow 传上游 pseudobulk 的 `outdir`，即
`diff/pseudobulk`）。列为 `gene, cell_type, comparison, baseMean, log2FC,
lfcSE, stat, pval, pval_adj, significant`；映射 `--gene-col gene --value-col
log2FC --p-col pval_adj`。`<tissue>` 目录与 `deg_<tissue>` 目录必须一致；
`all_pseudobulk_DEG.csv`（聚合表）不在扫描范围。单元发现由
`bin/discover_units.py` 在 checkpoint 规则里完成，路径知识只写在该引擎与
`.smk` 的 input 函数两处（同一契约）。

### 1.2 输出

`<outdir>/function/go_kegg/{counter}/{mode}/{tissue}/{contrast}/{cell_type}/`，
每单元 8 个文件（与 go-kegg.r 产物同构，零富集时也全部产出）：

`go_up.csv` `go_down.csv` `kegg_up.csv` `kegg_down.csv`
`go_back_to_back.<fmt>` `kegg_back_to_back.<fmt>` `up_genes.txt` `down_genes.txt`

图按 `formats` 每种格式一份（默认 png）。单元清单
`<outdir>/function/go_kegg/units.txt` 与完成哨兵
`<outdir>/function/go_kegg/enrichment.done`（node.py 注册的面板产物）落在
`go_kegg/` 根目录。

### 1.3 参数（Params.function.pseudobulk_go）

| 键 | 说明 | 默认 |
|---|---|---|
| `enabled` | node.py 读取；true 才注册 `enrichment.done` 进 DAG | 必填 |
| `root` | pseudobulk 根目录；空则用 `indir` | null |
| `modes` | 参与富集的模式（all/gene/te） | `["all","gene","te"]` |
| `species` | human / mouse / macaque（必填，不猜；node.py 生成期校验） | macaque |
| `lfc_cut` / `p_cut` | 上下调基因阈值 | 1 / 0.05 |
| `top` | 图中每方向条目数 | 10 |
| `formats` | 图格式，多种同时输出（合法值 png pdf svg ps eps tif tiff jpg jpeg；go-kegg.r 为 R 脚本，ggsave 设备集，pgf/raw/rgba 传入即报错）。旧单数键 `format` 仍受理 | `["png"]` |

### 1.4 直接调用引擎

```bash
Rscript modules/function/bin/go-kegg.r \
  -i <contrast>_<cell_type>_DEG.csv -o <unit_outdir> -s macaque \
  --gene-col gene --value-col log2FC --p-col pval_adj \
  --lfc-cut 1 --p-cut 0.05 --top 10 --format png
```

单元发现引擎：

```bash
python modules/function/pseudobulk_go/bin/discover_units.py \
  --root <pseudobulk_root> --modes '["all","gene","te"]' \
  --mode-dirs '{"all":"deg_pseudobulk","gene":"deg_pseudobulk_gene","te":"deg_pseudobulk_te"}' \
  --out units.txt
```

### 1.5 规则与运行

规则：`function_pseudobulk_go_index`（checkpoint，发现单元）→
`function_pseudobulk_go`（每单元富集）→ `function_pseudobulk_go_done`
（哨兵）→ `function_pseudobulk_go_result`（面板聚合点，仅 input）。

直接 snakemake（脚本化调用推荐，配置不经 run_module 物化）：

```bash
snakemake -s workflow/Omics/modules/function/pseudobulk_go/pseudobulk_go.smk \
  --configfile /path/to/pseudobulk_go.config.json \
  --cores 24 --resources mem_mb=900000 --sdm=apptainer \
  function_pseudobulk_go_result
```

run_module（需要显式 target；`--config` 会被同目录模板物化覆盖，参数请用
点号覆盖或直接改模板副本）：

```bash
python workflow/Omics/run_module.py \
  --snakefile modules/function/pseudobulk_go/pseudobulk_go.smk \
  --Params.function.pseudobulk_go.root=/path/to/pseudobulk_root \
  --Params.function.pseudobulk_go.species=macaque \
  -o /path/to/out -t 8 --sdm apptainer \
  --target function_pseudobulk_go_result
```

真实运行配置示例（luanchao 项目，macaque；输出落在
`<outdir>/function/go_kegg/`）：

```json
{
  "indir": "/data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/pseudobulk",
  "outdir": "/data/16T/luosg/output/luanchao_mulatta/scRNAseq",
  "logdir": "/data/16T/luosg/output/luanchao_mulatta/scRNAseq/log",
  "ROOT_DIR": "/data/16T/luosg/workflow/Omics",
  "env": {
    "env_dir": "/data/16T/Database/env",
    "function": "/data/16T/Database/env/function/function.sif"
  },
  "Params": {
    "function": {
      "pseudobulk_go": {
        "modes": ["all", "gene", "te"],
        "species": "macaque"
      }
    }
  }
}
```

注意：`--config` 指向的路径会被 run_module 的模板物化覆盖——配置放独立文件，
勿与模板同路径。

## 2. 结果解读

- **`{go,kegg}_{up,down}.csv`**：clusterProfiler 富集结果表，列 `ID`
  （GO/KEGG 条目号）、`Description`（条目名）、`GeneRatio`/`BgRatio`
  （命中基因/背景基因比例）、`pvalue`、`p.adjust`、`qvalue`、`geneID`
  （命中基因，`/` 分隔）、`Count`（命中数）。up/down 两个方向分别成表；
  无富集条目时表只有表头。
- **`up_genes.txt` / `down_genes.txt`**：参与富集的上/下调基因名（去重），
  即 §3 判定规则筛出的输入基因。
- **图像结构与元件含义（`*_back_to_back.<fmt>`）**：双向条形图，纵向排列
  富集条目（自上而下按 pvalue 升序取 `top` 条），横轴为 `-log10(pvalue)`；
  上调方向的条形画在轴的一侧、下调方向画在另一侧（背靠背），条目名经
  清理与换行后标注在纵轴。条形越长表示富集越显著。两侧均无条目时图中央
  显示 `No enriched terms` 占位文字。看图要点：两侧条目不对称属正常
  （up/down 基因数不同）；te 模式常见占位图（TE 名无法映射基因库）。

## 3. 算法原理

1. **上/下调基因划分**：对 DEG 表按阈值过滤——up：`log2FC >= lfc_cut` 且
   `pval_adj <= p_cut`；down：`log2FC <= -lfc_cut` 且 `pval_adj <= p_cut`；
   `pval_adj` 为 NA 的行剔除，基因名去重。
2. **GO/KEGG 富集**：上、下调基因分开做过表征分析（clusterProfiler
   `enrichGO`/`enrichKEGG` 超几何检验），物种注释库按 `species` 选择
   （org.Hs.eg.db / org.Mm.eg.db / org.Mmu.eg.db）；多重检验校正产生
   `p.adjust`/`qvalue`。
3. **空结果语义**：基因列表为空或无富集条目时按契约仍产出全部 8 个文件
   （空表 + 占位图），te 模式的 TE 名大多无法映射基因库，空表属预期结果。

## 4. 参考文献

- Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing
  biological themes among gene clusters. OMICS. 2012;16(5):284-287.
  doi:10.1089/omi.2011.0118
- clusterProfiler / org.Hs.eg.db / org.Mm.eg.db / org.Mmu.eg.db:
  Bioconductor annotation packages (https://bioconductor.org/packages/clusterProfiler)
- 本模块引擎：`modules/function/bin/go-kegg.r`（脚本头注明功能与依赖）
