# scanpy auto 模式：marker_supplement LLM 自动产出 — 设计方案

状态：已实施（开放问题的决定见 §6 括号内"已定"）
日期：2026-11-01

## 1. 现状

- `marker_supplement` 是人工输入：`Params.scanpy.<counter>.auto.marker_supplement = {tissue: json 路径}`（四个 config 面：`config/scRNAseq.json` + schema、`modules/scanpy/scanpy.json` + schema），`scanpy.smk` 在 scanpy_auto / scanpy_reannotate 透传 `-X`。
- 契约：`{cell_type: {positive: [...], negative: [...]}}`，**display-only** — 只追加到 marker 面板（negative 行标题 `(exp. absent)`），不进 LLM 提示词、不进特异性打分（scrnaseq-automode-qc.md 红线："do not silently change annotation behavior with it"）。
- LLM 注释已输出两组相关证据但不进 supplement：per-cluster `canonical_markers`（阳性）、`evidence.ruling_out`（阴性，按竞争类型分组，仅对 panel 表低 % 基因断言 absence）、`evidence.programs`（组合程序）。
- 结果：近缘细胞类型的判别组合 marker 全靠人工写 JSON，面板默认没有判别信息。

## 2. 目标

auto 模式注释收敛后，LLM 自动产出 per-cell-type 判别 marker 组合（阳性 + 关键阴性），持久化并自动进入 marker 面板；人工 JSON 保留为覆盖/补充输入；display-only 语义不变。

## 3. 产出机制选项

**A. 机械派生（0 额外 LLM 调用）**：从最终 annotation 聚合 — positive = `canonical_markers` ∪ `programs` 基因，negative = `ruling_out` 基因，按 cell_type 合并。免费、确定性强；但只是证据重排，不是"判别组合"精选，positive 易过泛，需要额外筛选规则。

**B. 专门 LLM 调用（推荐）**：注释收敛、`_auto_normalize_annotation_names` 之后，把最终 annotation 表 + 证据块交给 LLM，一次调用生成判别组合 `{cell_type: {positive, negative}}`（只输出有判别价值的类型），`var_names` 校验后原子写盘。与 `_canonicalize_cell_type_wording` 同模式（一次额外调用，失败保留空 supplement 不阻塞）。

**C. A+B**：机械派生兜底，LLM 精选覆盖。

## 4. 改动面（按 B 方案）

| # | 位置 | 改动 |
| --- | --- | --- |
| 1 | `modules/scanpy/bin/scRNAseq.py` | 新增 `_generate_marker_supplement()`：提示词 + JSON 解析 + `var_names` 校验（缺失基因 warning 跳过，同 `_load_marker_supplement` 规则），原子写 `<output>_auto_reports/table/marker_supplement.auto.json`；挂在 mode_auto 循环收敛后、`_auto_finalize_outputs` 之前（消费最终 annotations，而非中间迭代） |
| 2 | 同上（绘图路径） | `plot_markers` 的 `supplement_map` 改为合并来源：`marker_supplement.auto.json` + 人工 `-X`；mode_reannotate 自动发现同目录 `marker_supplement.auto.json` 合并 |
| 3 | 同上（CLI） | auto 子命令加 `--no-auto-marker-supplement` 关闭开关 |
| 4 | `modules/scanpy/scanpy.smk` | scanpy_auto 透传新开关（config bool） |
| 5 | 四个 config 面 | `Params.scanpy.<counter>.auto.marker_supplement_llm`（bool） |
| 6 | 文档 | ANNOTATION.md §2.3 补自动产出说明；scrnaseq-automode-qc.md 的 supplement 段同步（保持 display-only 红线措辞） |
| 7 | 测试 | 解析/合并逻辑用假 LLM 响应做单测；`run.py --test scRNAseq` 回归门 |

## 5. Non-goals

- 不把 supplement（无论人工还是 LLM 版）喂回提示词或特异性打分 — display-only 红线不动。
- 不改 `evidence.ruling_out` → `plot_evidence` 差分图链路。
- 不回填已有 run 的 supplement（新产出只对新跑的 auto 生效；旧 run 可走 no-op override + reannotate 图刷新路径）。

## 6. 开放问题（已全部决定）

1. **产出机制**：A / B / C？（已定：B，专门 LLM 调用）
2. **默认开关**：auto 模式默认开还是默认关？（已定：默认开，`--no-auto-marker-supplement` / `marker_supplement_llm: false` 关闭）
3. **合并策略**：人工 JSON 与 LLM 自动版冲突时，按 cell_type key 整条覆盖，还是 gene 级并集？（已定：key 级覆盖，大小写不敏感）
4. **生成范围**：只对 confusable sibling types 生成，还是全部注释类型？（已定：LLM 自判，只输出有判别价值的类型，`{}` 合法）
5. **改名失配**：reannotate 覆写改名后 supplement 的 cell_type key 失配怎么办？（已定：维持大小写不敏感匹配 + warning，不迁移 key）
