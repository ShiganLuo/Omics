# openms/proteininference — 蛋白推断模块(OpenMS Epifany)

OpenMS 工具链的蛋白推断步:对上游 FDR 过滤后的 idXML 用 Epifany
(Bayesian 蛋白推断)聚合肽段证据到蛋白/蛋白组,输出蛋白级 idXML。
子模块共享父目录 conda 环境(`../openms.yaml`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 上游 psm_fdr 结果目录(默认 `data/psm_fdr`);input 函数拼 `{indir}/{sample_id}/{sample_id}_filtered.idXML` |
| `samples` | 样本 ID 列表 |

输入 idXML 必须含肽段命中(`<PeptideHit`);零 PSM 存活时规则提前报错
(避免 Epifany 在空输入上段错误),指引先排查搜库/FDR 步。

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_protein.idXML   # 蛋白推断后 idXML
```

收尾规则 `protein_inference_result` 聚合全部样本输出。

### 1.3 参数

config 位于 `Params.protein_inference` 下,与 Epifany CLI 对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `method` | — | `epifany` | 推断方法;本 Epifany 构建仅支持 Bayesian `epifany`,其他取值直接报错 |
| `greedy_group_resolution` | `-greedy_group_resolution` | `none` | 蛋白组贪婪解组策略(TOPP 级参数) |
| `Procedure.epifany` | — | `Epifany` | Epifany 可执行路径 |
| mem_mb | — | 32000 | `rule_res("protein_inference")` 配置 |

直接调用示例(脱离流程):

```bash
Epifany -in S1_filtered.idXML -out S1_protein.idXML \
  -greedy_group_resolution none
```

## 2. 算法原理与生物学意义

蛋白质组搜库鉴定的是肽段;一个肽段可能对应多个共享同序列的蛋白
(蛋白组,protein group),蛋白推断把肽段证据按简约原则/概率模型归属
到蛋白。Epifany 用贝叶斯网络对 PSM 概率做蛋白层面的后验推断,并可
用 `-greedy_group_resolution` 对蛋白组做贪婪拆分,输出带推断结果的
idXML(`greedy_group_resolution=none` 保留蛋白组不拆)。生物学意义:
蛋白推断决定"鉴定到哪些蛋白"这一最终结论的可靠性,是定量(如
MSstats)前的必经步骤;结果读法:下游 quantification/报告读
`*_protein.idXML` 中的蛋白得分与蛋白组结构。

## 3. 参考文献

1. Pfeuffer J, Sachsenberg T, Dijkstra TMH, Serang O, Reinert K,
   Kohlbacher O. *EPIFANY: A Method for Efficient High-Confidence Protein
   Inference*. Journal of Proteome Research 19(3):1060-1072, 2020.
   DOI: 10.1021/acs.jproteome.9b00566
2. OpenMS — https://openms.de/;Röst HL, Schmitt U, Aebersold R,
   Malmström L. *OpenMS: a flexible open-source software platform for
   mass spectrometry data analysis*. Nature Methods 13(9):741-748, 2016.
   DOI: 10.1038/nmeth.3959
