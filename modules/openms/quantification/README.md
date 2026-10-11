# quantification（蛋白定量：LFQ / TMT / DIA）

从蛋白推断结果（`6_protein_inference/{sample}/{sample}_protein.idXML`）做蛋白定量，
按 `quantification_method` 分三条路线。方法选择由 subworkflow（node.py 注册输出）决定，
模块内定义全部规则。

## 规则

| 规则 | 输出（`<outdir>/`） |
|---|---|
| quantification_lfq | `lfq_quantification.mzTab` + `lfq_msstats.csv`（+ `lfq_quantification.triqler.tsv`，`Params.add_triqler_output` 开启时） |
| quantification_tmt_reporter | `{sample}/{sample}.reporter.consensusXML` / `{sample}.mapped.consensusXML`（每样本） |
| quantification_tmt | `tmt_quantification.mzTab` + `tmt_quantification.consensusXML` + `tmt_quantification.proteins.csv` |
| quantification_dia | `dia_quantification.mzTab` + `dia_msstats.csv`（+ `dia_quantification.triqler.tsv`） |
| quantification_result | 聚合目标（按 quantification_method 取对应 mzTab） |

## 方法原理

- **lfq/dia**：ProteomicsLFQ 一次产出 mzTab 和 MSstats 输入 CSV（`-out_msstats`）。
  蛋白汇总参数映射到 `-ProteinQuantification:*`（top:N / top:aggregate /
  best_charge_and_fraction）；`Params.protein_fdr` 映射到 `-proteinFDR` /
  `-picked_proteinFDR true`；`Params.quantify_decoys` 映射到
  `-PeptideQuantification:quantify_decoys`。
- **tmt**：报告离子定量链 IsobaricAnalyzer（`-type <tmt*>`，每样本提取 16/18 通道
  报告离子）→ IDMapper（把肽段鉴定映射到报告离子 consensus 特征）→ FileMerger
  （合并分级）→ ProteinQuantifier（mzTab 报告）。MSstats 级输入由 msstats 模块的
  MSstatsConverter（`-method ISO`）从合并后的 consensusXML 生成。

## config

```
expdesign                  OpenMS 实验设计 TSV（Params.workflow.expdesign）。LFQ 可选，
                           TMT 的 MSstats 转换必填
Params.protein_quant       quantification_method / protein_inference(aggregation|bayesian) /
                           top / average / best_charge_and_fraction
Params.tmt                 type(tmt6/10/11/16/18plex) / reporter_mass_shift /
                           isotope_correction / normalization / reference_channel /
                           correction_matrix / channel_descriptions {"126": "S1", ...}
Params.add_triqler_output  额外输出 Triqler 输入 TSV
Params.quantify_decoys     是否定量 decoy 肽段
Params.protein_fdr         fdr / picked（ProteomicsLFQ 蛋白水平 FDR）
```

旧键名兼容（自动映射）：`protein_quant.method`→quantification_method、
`best_charge_and_score`→best_charge_and_fraction、`tmt.channel_mass_tolerance`→
reporter_mass_shift、`tmt.channel_annotation`（"126:S1,127N:S2"）→channel_descriptions。

## 一键运行（经 subworkflow）

```
python run.py -m meta.tsv -w QuantMS -o <outdir> --sdm apptainer \
    --Params.workflow.quantification_method lfq
```

TMT（需提供实验设计）：

```
python run.py -m meta.tsv -w QuantMS -o <outdir> --sdm apptainer \
    --Params.workflow.quantification_method tmt \
    --Params.workflow.expdesign /path/to/expdesign.tsv \
    --Params.tmt.type tmt16plex
```

## 注意

- `quantification_method='tmt'` 时 searchengine 会先校验谱图含 TMT/iTRAQ 报告离子
  （`check_reporter_ions`），未标记数据会在搜库前明确报错，避免全噪声结果。
- MSstats 包不在 openms.sif 中；`8_msstats/msstats_input.csv` 是 MSstats R 包的
  ready-to-use 输入表，统计分析需另行安装 MSstats 后运行。
