# pseudobulk_go — scRNAseq pseudobulk 结果的 GO/KEGG 富集子模块

`modules/function` 的布局 C 子工具。扫描 scRNAseq pseudobulk 输出树，对每个
`(mode, tissue, contrast, cell_type)` 的 DEG 表做 up/down 分开的 GO/KEGG 富集，
引擎复用 `modules/function/bin/go-kegg.r`（clusterProfiler + org.Hs/Mm/Mmu.eg.db，
human/mouse/macaque）。te 模式同样受理（TE 名称无法映射基因库时产空表，属正常结果）。

## 输入结构（pseudobulk 输出契约）

```
<root>/<tissue>/pseudobulk/<counter>/deg_pseudobulk{,_gene,_te}/
    deg_<tissue>/<contrast>/<cell_type>/<contrast>_<cell_type>_DEG.csv
```

列为 `gene, cell_type, comparison, baseMean, log2FC, lfcSE, stat, pval,
pval_adj, significant`；映射 `--gene-col gene --value-col log2FC
--p-col pval_adj`。`<tissue>` 目录与 `deg_<tissue>` 目录必须一致；
`all_pseudobulk_DEG.csv`（聚合表）不在扫描范围。

## 输出

`<outdir>/function/go_kegg/{counter}/{mode}/{tissue}/{contrast}/{cell_type}/`，
每单元 8 个文件（与 go-kegg.r 产物同构，零富集时也全部产出）：

`go_up.csv` `go_down.csv` `kegg_up.csv` `kegg_down.csv`
`go_back_to_back.<format>` `kegg_back_to_back.<format>` `up_genes.txt` `down_genes.txt`

## 配置（Params.function.pseudobulk_go）

| 键 | 说明 | 默认 |
|---|---|---|
| `root` | pseudobulk 根目录；空则用 `indir` | null |
| `modes` | 参与富集的模式（all/gene/te） | `["all","gene","te"]` |
| `species` | human / mouse / macaque（必填，不猜） | macaque |
| `lfc_cut` / `p_cut` | 上下调基因阈值 | 1 / 0.05 |
| `top` | 图中每方向条目数 | 10 |
| `format` | 图格式 png/pdf/svg | png |

## 运行

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
  "indir": "",
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
        "root": "/data/16T/luosg/output/luanchao_mulatta/scRNAseq/common/5_combine_h5ad",
        "modes": ["all", "gene", "te"],
        "species": "macaque"
      }
    }
  }
}
```

注意：`--config` 指向的路径会被 run_module 的模板物化覆盖——配置放独立文件，
勿与模板同路径。
