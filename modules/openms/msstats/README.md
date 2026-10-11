# msstats（MSstats 输入表准备）

把定量结果转成 MSstats R 包的输入 CSV，输出 `<outdir>/msstats_input.csv`。
（该文件是 MSstats 的输入而非统计结果；R MSstats 包不在 openms.sif 中，
统计分析 dataProcess/groupComparison 需在装有 MSstats 的环境中运行。）

## 规则

| 规则 | 输出 |
|---|---|
| msstats | `msstats_input.csv` |
| msstats_result | 聚合目标 |

## 方法原理

- **lfq/dia**：ProteomicsLFQ 的 `-out_msstats` 已经是 MSstats 输入格式，
  本步骤只把它暂存到统一的最终路径。
- **tmt**：MSstatsConverter `-method ISO` 从合并的报告离子 consensusXML
  （`7_quantification/tmt_quantification.consensusXML`）+ 实验设计 TSV 生成
  MSstatsTMT 输入表。`Params.workflow.expdesign` 必填，缺失时报错并指明 JSON 键。

## config

```
Params.workflow.expdesign  OpenMS 实验设计 TSV（TMT 必填）
quantification_method      lfq / tmt / dia（决定消费哪种定量产物）
```

## 运行（经 subworkflow）

```
python run.py -m meta.tsv -w QuantMS -o <outdir> --sdm apptainer \
    --Params.workflow.quantification_method tmt \
    --Params.workflow.expdesign /path/to/expdesign.tsv
```

输出：`<outdir>/common/8_msstats/msstats_input.csv`。
