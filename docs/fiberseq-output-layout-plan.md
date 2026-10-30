# Fiberseq 输出目录重构方案 v4（已执行）

状态：**已执行完毕（2026-10-30）**。代码、存量迁移、DAG 验证全部完成。
v4：输出目录只保留流程产物。v3：`results/` → `fiberseq_analysis/`。
v2：放弃「variation 等收编」「work/ 工作区」。
实测对象：`output/Totipotency/Fiberseq/`（共 1.1T）。

## 一、目标布局（纯流程产物，别无它物）

```
Fiberseq/                       # 流程输出目录，只含流程产物
├── raw.json                    # 运行配置快照
├── common/                     # 中间大文件，编号步骤，每文件唯一路径
│   ├── 1_raw_data/<s>/         # 原始 HiFi BAM 符号链接
│   ├── 2_fiberseq_bam/<s>/     # ft predict-m6a / add-nuc / fire
│   ├── 3_align_bam/<s>/        # pbmm2 sorted BAM —— 全流程唯一排序 BAM
│   ├── 4_raw_results/<s>/      # 元素 BED 唯一权威路径：m6a/nuc/msp/cpg.bed.gz、
│   │                           # fire_peaks.bed、fire_fdr_table.tsv、qc.tsv
│   │                           # （hap{1,2} 四件套 BED 也在这里）
│   ├── 5_phased_bam/<s>/       # hiphase phased BAM + hap1/hap2 BAM —— 唯一 hap BAM
│   └── 6_coverage_bw/<s>/      # deeptools bigWig
├── fiberseq_analysis/          # fiberseq 分析家族交付物（表/图）
│   ├── annotate/<s>/  motif/<s>/  footprint/<s>/  nucpos/<s>/
│   ├── quant/<s>/ + quant/matrix/
│   ├── track/<s>/  coactuation/<s>/  haplotype/<s>/  censat/<s>/
│   └── diff/<contrast>/
├── variation/  eccdna/  function/  ldsc/  genome/   # 顶层，不动
├── log/                        # sample/ group/ + cmd/（命令回执）
├── tmp/                        # 纯流程 tmp_dir，跑完即空
└── .snakemake/                 # 原位（snakemake 以 outdir 为 cwd，不动）

非流程产物的去处（全部移出本目录树）：
├── 手工分析（原 analysis/，4G）   → ../Fiberseq_manual/   # 与 Fiberseq 平级
└── 回归切片（tmp/realslice）      → workflow/Omics/test/data/fiberseq_realslice/
    测试残留（replay_check 等 77G）→ 直接删除
```

命名理由：顶层 variation/eccdna 是领域命名，`results` 过泛且与它们语义
打架；目录内 10 个子目录全部 `fiberseq_` 前缀属冗余。`fiberseq_analysis/`
与顶层领域目录平级一致，与 docs/plan、模块 README「下游分析家族」叫法
统一；改名后路径反而更短（`fiberseq_analysis/quant/DMSO/`）。

## 二、三个真问题与对应修复（流程可控的部分）

### 问题 1：同内容双份实体 ≈ 270G（md5 已验证）

| 重复对 | 修复 |
|---|---|
| hap BAM：`results/<s>/` 与 `common/4_raw_results/<s>/` 各一份（≈41G） | 根因是 `node.py runFiberseq` 注册 `results/<s>/<s>.hap*.bam`，而 `Fiberseq.smk` 的 `fibertools_hap_config.outdir` 指向 `common/4_raw_results`，两套路径各跑一遍（10-25 与 10-29 各生成一份）。改法：`fibertools_hap_config.outdir` → `common/5_phased_bam`（hap BAM 与 phased BAM 同居），`node.py` 目标同步改；`results/<s>/` 目录整个消失 |
| `fire.sorted.bam`：`3_align_bam` 与 `4_raw_results` 各一份（≈124G） | `ft_call_peaks` 在自己 outdir 又生成一份排序 BAM，无下游引用。改法：该规则直接用 `3_align_bam` 的输入 BAM，不再复制 |
| `fire.bed` + `fire.bed.orig` 双份未压缩（≈104G+46G） | `.orig` 是纯备份，删；`fire.bed` 是否转 bgzip 见开放问题 1 |

### 问题 2：元素文件无唯一权威路径（隐性 bug，已实际发生过）

`ft extract` 写 `common/4_raw_results/<s>/`，而 10 个 fiberseq_* 子模块的
`indir` 指向 `results/`（即读 `results/<s>/<s>.fire_peaks.bed`）。10-15 的
quant 日志显示当时读的是 `results/DMSO/DMSO.msp.bed.gz`，该文件如今只存在于
`common/4_raw_results/`——中途靠手工搬文件才对上。改法：
`Fiberseq.smk` 给分析模块族传 `indir = common/4_raw_results`；
footprint/diff/haplotype 的跨模块输入（motif_dir、quant_dir、annot_dir）
已有独立 config 键，显式指向 `fiberseq_analysis/*`，不随 indir 漂移。

### 问题 3：散落物
- 命令回执 `*_run_<s>_<时间戳>.sh` 混在结果目录（同一步多份）→ 统一写
  `{logdir}/cmd/{sample_id}/`（规则自身写文件，完全可控）。
- `.tmp` 中间残留（如 `motif_footprint_scores.tsv.tmp`）→ 成功后清理。
- 容器缓存 `.vep/.cache/.config/.keras/.parallel/.RepeatMaskerCache` 落在
  输出根：按规则注入环境变量（MPLCONFIGDIR、VEP_CACHE_DIR、
  XDG_CACHE_HOME，工具认哪个设哪个，指向 `tmp/` 下子目录）——只改产生
  这些缓存的规则（fiberseq 画图脚本、sv 的 VEP/RepeatMasker）；存量目录
  迁移时一次清掉。做不到的不硬控，README 记一笔。

## 三、工程改动清单

| # | 文件 | 改动 |
|---|---|---|
| 1 | `node.py runFiberseq` | 分析目标注册 `res = {outdir}/fiberseq_analysis`，子目录去前缀（`annotate/`…）；hap 目标 → `common/5_phased_bam/<s>/`；元素 BED 目标核对为 `common/4_raw_results` |
| 2 | `subworkflow/Fiberseq.smk` | `analysis_outdir = {outdir}/fiberseq_analysis`；各分析模块 `outdir` 去前缀；`fibertools_hap_config.outdir` → `common/5_phased_bam`；模块族 `indir` → `common/4_raw_results`；motif_dir/quant_dir/annot_dir 显式指向 `fiberseq_analysis/*`；function 的 indir（读 diff 基因表）同步 |
| 3 | `modules/fiberseq/*/*.smk`（10 个） | 默认 `outdir` 改 `fiberseq_analysis/<name>`（无前缀）；回执入 `log/cmd/{sample}/`；`.tmp` 清理；画图规则加 MPLCONFIGDIR |
| 4 | `modules/fibertools/fibertools_analysis.smk` | `ft_call_peaks` 不再复制排序 BAM；shuffle 空床写 tmp_dir；回执入 `log/cmd/` |
| 5 | sv/decoil 相关规则 | VEP_CACHE_DIR、RepeatMasker 缓存指向 tmp/ 子目录 |
| 6 | `fiberseq.schema.json`、模块 README、`docs/fiberseq-downstream-files.md` | 路径约定、输出速查同步（`fiberseq_analysis/<模块>/<样本|对比>/`） |

不改：`.snakemake` 位置、tmp/ 名称、variation/eccdna/function/ldsc 位置、
run.py 全局行为。

## 四、存量迁移（一次性脚本：先输出 mv/rm 清单 → 人工过目 → 执行）

流程目录内：

| 动作 | 内容 | 释放 |
|---|---|---|
| rm | `results/<s>/<s>.hap*`（与 common 副本 md5 一致已验证） | ≈41G |
| mv | `common/4_raw_results/<s>/<s>.hap*` → `common/5_phased_bam/<s>/` | 0 |
| rm | `common/4_raw_results/<s>/<s>.fire.sorted.bam(.bai)`（无下游引用） | ≈124G |
| rm | `<s>.fire.bed.orig` × 2 | ≈104G |
| mv | `results/fiberseq_<name>` → `fiberseq_analysis/<name>`（10 个 + matrix 随 quant） | 0 |
| rm | `.vep .cache .config .keras .parallel .RepeatMaskerCache`（可再生缓存） | ≈0.9G |
| rm | 各结果目录散落回执 `.sh`、`.tmp`（复制到 log/cmd/ 留档后删） | 0 |

移出目录树（非流程产物）：

| 动作 | 内容 | 去处 |
|---|---|---|
| mv | `analysis/`（4G 手工分析） | `../Fiberseq_manual/`（与 Fiberseq 平级） |
| mv | `tmp/realslice`（回归切片） | `workflow/Omics/test/data/fiberseq_realslice/` |
| rm | `tmp/{replay_check, res_design_demo, rewrite_test, track_test, tmp_*, phased_reads_*.tsv}` | 删除（≈77G） |
| mv→rm | `tmp/__pycache__` 等杂项 | 删除 |

预计净释放 ≈340G。迁移后 `--rerun-triggers mtime` 或 touch 校正元数据
（新路径文件沿用旧 mtime，避免误判重算）。

## 五、验证

1. 全开关 `--dry-run`：DAG 无 MissingInput；无目标落在旧 `results/` 路径。
2. realslice 真实小切片回归（切片位于
   `workflow/Omics/test/data/fiberseq_realslice/`）：新路径跑
   extract→peaks→quant→motif，与旧结果 checksum/行数比对。
3. 迁移前后 `du` 对账；`fiberseq_analysis/` 总量 <10G 且不含 BAM/大 BED。

## 六、执行结果（2026-10-30）

开放问题 1/2/3 均经确认（bgzip 化、更名 4_elements、手工与切片外迁）。

- 代码：subworkflow/Fiberseq.smk、node.py、10 个 fiberseq 子模块 smk+json、
  fibertools_analysis.smk、quant/coactuation fire 输入改 .gz、
  fiberseq_track.py 成功转 bw 后删 bedGraph、motif/footprint 启动清 .tmp。
- 迁移：127 项动作（fiberseq_migration.log），目录 1.1T → 608G，磁盘
  52% → 49%（净释放约 400G）。tabix 抽查 fire.bed.gz 首尾记录正常。
- 验证：新 raw.json 126 个目标全部存在于新路径；真实目录 dry-run DAG
  干净，仅 8 个轻量任务（quant/track/coactuation/merge）因输入 mtime
  变化待重跑——这是 gz 化后的端到端验证点。
- 归档：43+ 份历史回执入 `log/cmd/`；`analysis/` → `../Fiberseq_manual/`
  （4G）；realslice → `workflow/Omics/test/data/fiberseq_realslice/`；
  77G 测试残留与 0.9G 缓存目录已删。
