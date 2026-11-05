# modules/sv — 分析方法文档模块

本模块收录 sv 模块各分析线的**方法文档**(与代码同版本管理, 供审核/复现)。
文档只记录方法: 输入契约、算法、判定规则、统计设计、参数依据、质控点、
方法学教训、局限。**不含具体数据结果**(结果在数据输出目录)与
**运行用法**(用法见代码 argparse 帮助与 sv.smk 规则调用)。

| 文档 | 覆盖分析 | 对应代码入口 |
|---|---|---|
| `MEI_METHODS.md` | MEI 管线方法(输入契约/五步流程/判定逻辑三层/局限) | `bin/mei_analysis.py`(**唯一入口**, 分析→图一条龙; 绘图经 `bin/utils/mei_plot.py` 自动串联) |
| `README_polyA_TSD.md` | poly(A)/TSD 指纹检测算法 + TPRT 生物学原理 | `bin/mei_analysis.py` 的 fingerprint 步骤 |
| `TRA_TE_ANALYSIS_METHODS.md` | TRA/TE-TE 重组与年轻插入验证(T0-T3 证据漏斗、统计设计、VAF 验证、SURVIVOR 行为定性) | `bin/utils/tra_te_analysis.py` + `bin/utils/te_age_analysis.py` |

约定:
- 结论标注证据级(T0-T3 等), 不确定的写"候选/假设", 不下硬结论;
- 参数给出默认值与依据;
- 历次方法学错误的修正记录保留在文档内("方法学教训"节), 供审核追溯。

运行环境见 `../sv.yaml` / `../sv.def`(sv.sif 容器: RepeatMasker/blastn/
pysam/pandas/matplotlib); samtools 用 `/data/16T/Database/env/samtools/samtools.sif`。
