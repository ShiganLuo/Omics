# openms/psmfdr — PSM 级 FDR 控制子模块(OpenMS FalseDiscoveryRate)

对打分后的 PSM(psmrescoring 产出的 idXML)做 target-decoy FDR 控制,
输出过滤后的 idXML。conda 环境共享父目录 `../openms.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config,默认 `data/psm_rescoring`)+ 样本 ID 布局(形态 A) | `{indir}/{sample_id}/{sample_id}_scored.idXML`:psmrescoring(Percolator 打分)产出的 idXML |
| `samples` | 样本 ID 列表(供 subworkflow 注册驱动 DAG) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_filtered.idXML   # FDR 过滤后的鉴定结果
{logdir}/{sample_id}/psm_fdr.log
```

收尾规则 `psm_fdr_result` 聚合过滤后 idXML。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.openms` | 可执行文件本身 | `FalseDiscoveryRate` | OpenMS FDR 工具可执行路径 |
| `Params.psm_fdr_control.fdr` | `-FDR:PSM` | 0.01 | PSM 级 FDR 阈值 |
| `Params.psm_fdr_control.method` | — | `BH` | `BH`/`qvalue`:Benjamini-Hochberg q 值(默认);`strict`:追加 `-algorithm:no_qvalues`,朴素 FDR 不做 q 值平滑;其余值报错 |
| `indir` / `outdir` / `logdir` | — | `data/psm_rescoring` / `output` / `logs` | 输入/输出/日志根目录 |

资源默认:`mem_mb = 16000`(经 `rule_res("psm_fdr")` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
FalseDiscoveryRate \
  -in s1_scored.idXML -out s1_filtered.idXML \
  -FDR:PSM 0.01 -protein false          # method=strict 时追加 -algorithm:no_qvalues
```

## 2. 算法原理与生物学意义

蛋白组鉴定中 target-decoy 策略:把 decoy 序列混入搜索,以 decoy 命中率
估计假阳性率。OpenMS `FalseDiscoveryRate` 从 Percolator 分数重算
target-decoy q 值并覆盖原 q 值注释——它是本链路唯一的 FDR 权威(过滤与
报告的 q 值同源,避免两套估计器打架);`-protein false` 只控制 PSM 级
FDR。`method=BH` 用 Benjamini-Hochberg 平滑出 q 值,`strict` 用朴素 FDR。
生物学意义:FDR 控制决定鉴定结果的可信度,是定量与差异分析的前提;1%
PSM FDR 是领域常规。结果阅读:过滤后 idXML 中每 PSM 带 q 值注释,
`protein=false` 下蛋白推断由下游(psmrescoring/proteininference 链)负责。

## 3. 参考文献

1. OpenMS —— https://www.openms.org/(FalseDiscoveryRate 工具文档
   https://openms.readthedocs.io/)
2. Röst HL, et al. *OpenMS: a flexible open-source software platform for
   mass spectrometry data analysis*. Nature Methods 2016, 13(9):741–748.
   DOI: 10.1038/nmeth.3959
