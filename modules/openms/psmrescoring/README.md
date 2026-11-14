# openms/psmrescoring — PSM 特征提取与 Percolator 重打分子模块

C 形态子模块(单规则 `psm_rescoring` + result 锚点):对搜库产出的 idXML
先做 PSMFeatureExtractor 特征补全,再用 OpenMS PercolatorAdapter 重打分
(半监督 SVM,控制 FDR)。conda 环境共享父目录 `../openms.yaml`(见
`modules/openms/`)。

## 1. 模块输入输出与参数

### 1.1 输入

按 `Params.search_engine` 选择搜库引擎命名的 idXML(形态 A):

| 输入 | 条件 |
|---|---|
| `{indir}/{sample_id}/{sample_id}_comet.idXML` | `search_engine == "comet"` |
| `{indir}/{sample_id}/{sample_id}_msgf.idXML` | `"msgf"` |
| `{indir}/{sample_id}/{sample_id}_sage.idXML` | `"sage"` |
| `{indir}/{sample_id}/{sample_id}_msfragger.idXML` | `"msfragger"` |

`search_engine` 取其他值时抛 `ValueError`。

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_features.idXML     # 特征补全后的中间 idXML
{outdir}/{sample_id}/{sample_id}_scored.idXML       # Percolator 重打分结果
{outdir}/{sample_id}/psm_rescoring.sh               # 回执(两条命令)
logdir/{sample_id}/psm_rescoring.log                # 运行日志
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.search_engine` | — | — | 搜库引擎(comet/msgf/sage/msfragger),决定输入名 |
| `Procedure.percolator` | 可执行文件 | null(回退 `PercolatorAdapter`) | OpenMS PercolatorAdapter |
| `Procedure.psm_feature_extractor` | 可执行文件 | null(回退 `PSMFeatureExtractor`) | OpenMS PSMFeatureExtractor |
| `Params.psm_rescoring.percolator.train_FDR` | `-trainFDR` | 0.05 | 训练集 FDR 阈值 |
| `Params.psm_rescoring.percolator.test_FDR` | `-testFDR` | 0.05 | 测试集 FDR 阈值 |
| `Params.psm_rescoring.psm_feature_extractor.extra` | 附加参数串 | "" | 按空格拆分追加给 PSMFeatureExtractor |
| `samples` | — | [] | 样本列表(result 锚点 expand) |
| `indir` / `outdir` / `logdir` | — | `data/search_engine` / `output` / `logs` | 目录基路径 |

注意:json 中 `Params.psm_rescoring.percolator.feature`(`"top_psm"`)为
OpenMS 2.x 遗留键,3.x 的 PercolatorAdapter 已无 `-feature` 选项,.smk 不
读取。OpenMS 3.x 选项为 `-trainFDR/-testFDR`(无下划线)。

内存经 `rule_res("psm_rescoring")` 配置:mem_mb 默认 32000。

引擎命令(回执脚本原文):

```bash
PSMFeatureExtractor -in {s}_comet.idXML -out {s}_features.idXML [ extra... ]
PercolatorAdapter -in {s}_features.idXML -out {s}_scored.idXML -trainFDR 0.05 -testFDR 0.05
```

## 2. 算法原理与生物学意义

纯工具封装:肽谱匹配(PSM)重打分分两步——PSMFeatureExtractor 把搜库
引擎特有的打分特征(如 Comet E-value/Xcorr、MS-GF:SpecEValue)补进
idXML,Percolator 再以半监督 SVM 综合特征重打分,并按 `-trainFDR`/
`-testFDR` 控制错误发现率。生物学意义:搜库引擎原始分数校准差、可比性
弱;统一 Percolator 重打分后 PSM 的 q 值可跨引擎比较,是控制鉴定假阳性、
产出可信蛋白/肽段列表的关键步骤。

## 3. 参考文献

1. OpenMS(PercolatorAdapter / PSMFeatureExtractor)—
   https://www.openms.org/ ;https://openms.readthedocs.io/
2. Percolator — Käll L et al. *Semi-supervised learning for peptide
   identification from shotgun proteomics datasets*. Nat Methods 4(11):923–925,
   2007. DOI: 10.1038/nmeth1113
