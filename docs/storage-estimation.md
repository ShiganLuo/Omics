# 运行数据大小估算（storage estimation）

[← 返回 README](../README.md)

**一句话**：跑之前用 `estimate_storage.py estimate` 估算一次运行会占用多少磁盘（最终保留 / 运行期暂用 / 峰值三个数），跑完用 `measure`/`calibrate` 实测校准，越跑越准。

工具：[`src/common/util/estimate_storage.py`](../src/common/util/estimate_storage.py)（stdlib-only CLI，与 `calibrate_resources.py` 同风格：append-only 历史、不猜、区间输出、可审计）。

---

## 1. 背景

PacVar 等长读长工作流单次运行的输入可达数百 GB，BAM 中间产物成倍放大。在提交运行前需要回答：

- 跑完后这块盘要留多少（最终保留）？
- 运行过程中最多会占到多少（峰值，须按此预留）？
- 哪些大文件是中间产物、会自动删（暂用）？

手工估算（输入字节 × 经验倍数）容易漏：哪些输出被 `temp()` 标记、config 开关关掉了哪些步骤、峰值出现在哪个阶段。脚本把这套知识代码化，并用真实运行持续校准。

## 2. 用法

```bash
cd workflow/Omics

# ① 跑之前估算（核心功能）
python src/common/util/estimate_storage.py estimate \
    --meta data/meta_PacVar_iteration_totipotency.tsv \
    --config config/PacVar.json --workflow PacVar --json report.json

# ② 跑完（或跑到一半）实测 —— 任意工作流，不需要产物清单
python src/common/util/estimate_storage.py measure --root OUTDIR

# ③ 用真实运行校准经验比例（append-only）
python src/common/util/estimate_storage.py calibrate \
    --root OUTDIR --meta data/meta_xxx.tsv --workflow PacVar
```

`estimate` 读取校准历史（默认脚本同目录 `storage_ratio_history.tsv`，`--calibration` 可改），
同类产物有 ≥3 个实测点后自动用实测分位数覆盖经验区间。

## 3. 估算原理（四层）

### 3.1 输入基准：真实 stat，不猜

meta 里每个输入文件（`fastq_1/2`、`bam`/`pbi`、`ms_file` 列）真实 `stat` 取字节数；
文件不存在直接报错（`--allow-missing` 可放行，但会低估）。样本大小不同（如 4 个
HiFi BAM 为 108–140 GB 不等）就逐样本算再求和，不用均值 × n。

输入数据通过软链接进 `common/1_raw_data/`（MetaUtil `_link_file`），不占额外空间，
估算不计入。

### 3.2 产物清单（manifest）：每类产物 = 基准 × 比例

每个工作流一份产物类清单（`MANIFESTS`，PacVar 已内置），每个类声明：

| 字段 | 含义 |
|---|---|
| `scope` | `sample`（每样本一份 × n）/ `run`（全局一份） |
| `basis` | 大小基准：`input_bam` / `input_fastq` / `input_total` / `genome_fasta` / `fixed` |
| `ratio` | (low, mid, high) 三档：非 fixed 为「输出字节/基准字节」比例；fixed 为每单位绝对字节 |
| `gate` | config 开关条件（`Params.<step>.enabled` 等），不满足则该类不产生 |
| `temp_when` | 满足时该类是 temp（跑完 Snakemake 自动删）|
| `born`/`die` | stage 时间线：产生阶段 / 最后一个消费者所在阶段（temp 在该 stage 仍占用）|
| `patterns` | measure/calibrate 归类用的路径 glob |

比例系数的语义与取值依据（详见第 6 节精度边界）：

- BAM 类 ≈ 1.0× 输入 bam：pbmm2/GATK 只重新压缩同一批读段，±5%；
- bigWig / VCF / 组装目录按同类项目经验给宽区间；
- fixed 类（QC 表、图、报告）给绝对字节三档。

### 3.3 kept / temp 语义：与模块代码一致

产物是否保留不是拍脑袋，而是逐一对齐模块的输出声明：

- `maybe_temp(..., "gatk_prepare")` 包裹的产物（如 `addReadsGroup.bam`）在
  `Params.gatk_prepare.keep_intermediates == False` 时标 `temp()`，消费者跑完即删；
- 裸 `output:`（如 pbmm2 的 `{s}.sorted.bam`、MarkDuplicates 的
  `{s}.sorted_markdup.bam`）**不受 keep_intermediates 控制，永久保留** ——
  「不保留中间文件」实际只删 1 份 BAM，最终留 2 份，这点容易误估；
- config 的 `Params.<step>.enabled` 开关按 node.py `runPacVar` 的 gate 语义过滤
  （`_enabled` 缺省 True，sawfish 缺省 False），关掉的步骤不产生输出也不占空间。

### 3.4 峰值：stage 时间线

每个类带 `born`/`die`，峰值 = 各 stage 在盘字节之和的最大值：

- kept 类（`die=None`）：`born` 之后一直占用；
- temp 类：占用到 `die`（含）——Snakemake 在最后一个消费者完成前不删，
  所以 temp 在其消费者所在 stage 仍计数。

以 PacVar 默认配置为例（4 样本，输入 490 GB）：

```
stage 2  pbmm2 sorted.bam                       1×
stage 3  + addReadsGroup.bam                    2×
stage 4  + sorted_markdup.bam（temp 此后回收）    3×  ← 峰值
stage 5  + SV VCF / bigWig（约 2× + 小文件）
```

## 4. 输出解读

以 `meta_PacVar_iteration_totipotency.tsv`（4 个 HiFi BAM，合计 489.9 GB）+
`config/PacVar.json`（snv/telomere/repeat/phase/sawfish 关，sv/decoil 开）为例：

```
============================================================
最终保留(跑完后长期占用): 984.4 GB  [931.8 GB ~ 1.1 TB]
运行期暂用(temp,之后自动删): 489.9 GB  [465.4 GB ~ 538.8 GB]
峰值(stage 4,须按此预留): 1.5 TB  [1.4 TB ~ 1.6 TB]
============================================================
```

- **最终保留** ≈ 2× 输入：sorted.bam + sorted_markdup.bam 两份永久 BAM + 小文件；
- **暂用** ≈ 1× 输入：addReadsGroup.bam（keep_intermediates=false，MarkDuplicates 后回收）；
- **峰值** ≈ 3× 输入：三份 BAM 在 stage 4 短暂共存，**按此预留磁盘**；
- 区间 [low ~ high] 是比例系数的不确定性（校准后收窄）；
- 明细表列出每个类的 kept/temp、校准点数（`校准 n=N` = 实测比例，
  `经验比例 n=0` = 尚无实测）以及被 config 关闭的类。

`--json` 输出机器可读报告（classes / final_bytes / transient_bytes / peak_bytes /
peak_stage / unmodeled），供上层调度或容量规划消费。

## 5. 校准闭环（越跑越准）

`calibrate` 把一次真实运行折算成比例点，追加进
`storage_ratio_history.tsv`（append-only，同 `calibrate_resources.py` 的
`resource_history.tsv` 风格）：

```
run_id  ts  workflow  class_id  scope  sample_id  basis_bytes  observed_bytes  ratio  source
```

- 比例类记 `ratio = 实测字节 / 基准字节`；fixed 类记每单位绝对字节；
- 样本级产物按路径分量精确归属到 sample_id（归属不了的丢弃，不猜）；
- `estimate` 读历史：同类 ≥3 个点取 p50 为中值、p10/p90 为区间，覆盖经验值；
- 跨机器/跨数据类型积累后，比例自然分化（如 Fiberseq BAM 与普通 HiFi BAM 的
  比对放大系数不同），无需改代码。

建议节奏：每个新工作流/新数据类型首跑后 calibrate 一次，之后 estimate 即为实测锚定。

## 6. 精度边界（显式不猜）

- **区间而非单点**：经验比例给 low/mid/high 三档；只有校准点足够时才收窄；
- **显式未建模**（报告尾行注明）：GATK/SnakeMake 的 java tmp-dir 运行期临时写入
  （`out/tmp`，结束后清理）、小文件多时的文件系统块开销；
- **无清单的工作流**：`measure`/`calibrate` 走文件后缀兜底归类（bam/vcf/bigwig/
  fasta/...），`estimate` 不凭空编造产物——宁可报「无该工作流的产物清单」；
- 经验比例的依据：BAM 类来自重压缩不改读段的确定性推理（±5%）；
  bigWig/VCF/组装类来自同类项目观测，属待校准近似（与资源平衡提案 §1
  的自查条款同一口径：一切魔法数须可被校准闭环重定）。

## 7. 扩展新工作流

往 `MANIFESTS` 加一份产物类清单（字段见 §3.2 与脚本头 ClassSpec 注释），
或 `--manifest extra.json` 合并（格式 `{"WorkflowName": [类 dict, ...]}`）。

清单内容的依据必须是模块代码（输出声明 + temp 标记 + node.py 的 outfile
gate），不按「一般流程都有的东西」臆造。新增后跑一次 `calibrate` 积累实测点。
