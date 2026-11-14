# diann — DIA 蛋白组搜谱与定量模块(DIA-NN)

DIA-NN 全队列一次搜谱(library-free 或谱库模式)+ 自研汇总脚本:
输出肽段/蛋白鉴定报告、统计、每 run QC 表与蛋白定量矩阵。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | mzML 目录(上游 raw2mzml 的 `outdir`):`{indir}/{sample}/{sample}.mzML`(全部样本作为一个整体输入 diann_run) |
| `samples` (config) | 样本 ID 列表 |
| `genome.protein_fasta` (config) | 蛋白参考 FASTA(`--fasta`);**必填**,缺失时 `ValueError` |
| `Params.diann.library` (config) | 可选谱库(.tsv/.speclib);非空走库模式 |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 输出/日志/项目根目录 |

设计要点(见规则 docstring):全队列单 job 联合搜谱——DIA-NN 跨 run
对齐并在队列层面控制 FDR,是官方推荐用法。

### 1.2 输出

```
{outdir}/diann.report.tsv               # 主报告(肽段/蛋白鉴定与定量)
{outdir}/diann.report.stats.tsv         # 搜谱统计
{outdir}/diann_qc.tsv                   # 每 run QC 汇总(diann_summary.py)
{outdir}/diann_protein_matrix.tsv       # 蛋白 × 样本定量矩阵
```

回执:`{outdir}/diann_run.sh`、`diann_summary.sh`;日志
`{logdir}/group/diann_run.log`、`diann_summary.log`。

### 1.3 参数

config 位于 `Params.diann`、`Procedure`、`genome` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.protein_fasta` | `--fasta` | — | 必填 |
| `Params.diann.qvalue` | `--qvalue` | 0.01 | 鉴定 q-value 阈值 |
| `Params.diann.matrix_qvalue` | `--matrix-qvalue` | 0.01 | 矩阵 q-value 阈值 |
| `Params.diann.mass_acc_ms1` | `--mass-acc-ms1` | 20 | MS1 质量精度(ppm) |
| `Params.diann.mass_acc` | `--mass-acc` | 20 | MS2 质量精度(ppm) |
| `Params.diann.min_fr_mz` / `max_fr_mz` | `--min-fr-mz` / `--max-fr-mz` | 200 / 1800 | 前体 m/z 范围 |
| `Params.diann.min_pept_len` / `max_pept_len` | `--min-pept-len` / `--max-pept-len` | 7 / 30 | 肽段长度范围 |
| `Params.diann.fixed_mod` | `--fixed-mod` | `Carbamidomethyl (C),57.0215` | 固定修饰 |
| `Params.diann.var_mods` | `--var-mods` | `Oxidation (M),15.9949` | 可变修饰 |
| `Params.diann.library` | `--lib` | null | 非空 → 库模式,不再加 `--gen-spec-lib`/`--predictor` |
| `Params.diann.gen_spec_lib` | `--gen-spec-lib` | true | 无库时生成谱库 |
| `Params.diann.predictor` | `--predictor` | true | 无库时用神经网络预测谱库(library-free) |
| `Params.diann.verbose` | `--verbose` | true | 详细日志 |
| `Params.diann.extra` | 原样附加 | "" | 透传任意 DIA-NN 参数 |
| `Procedure.diann` | 引擎可执行 | `diann` | — |
| (threads) | `--threads` | 16 | `rule_res("diann_run")` 可覆盖;示例 config 的 `Params.diann.threads`/`enabled` 不被 .smk 消费 |

自研汇总脚本参数(`diann_summary.py`):

| 参数 | 默认 | 说明 |
|---|---|---|
| `-i` / `--input` | 必填 | DIA-NN report.tsv |
| `-o` / `--outdir` | 必填 | 输出目录 |
| `-q` / `--qvalue` | 0.01 | 鉴定 q-value 阈值 |
| `-l` / `--log` | null | 日志文件 |
| `-n` / `--dry-run` | false | 只记录计划不执行 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
diann --out out/diann.report.tsv --fasta protein.fasta \
    --threads 16 --qvalue 0.01 --matrix-qvalue 0.01 \
    --mass-acc-ms1 20 --mass-acc 20 \
    --min-fr-mz 200 --max-fr-mz 1800 \
    --min-pept-len 7 --max-pept-len 30 \
    --fixed-mod "Carbamidomethyl (C),57.0215" \
    --var-mods "Oxidation (M),15.9949" \
    --gen-spec-lib --predictor --verbose \
    s1.mzML s2.mzML

python modules/diann/bin/diann_summary.py \
    -i out/diann.report.tsv -o out -l log/diann_summary.log
```

## 2. 算法原理与生物学意义

DIA(数据非依赖采集)把每个扫描窗口内的全部前体碎裂,蛋白组覆盖
远超 DDA;DIA-NN 用神经网络预测谱/保留时间做 library-free 搜谱,并以
「信号校准 + q-value 控制」在跨 run 对齐后整体控制 FDR(所以本模块
全队列一次跑而非逐样本)。关键阈值解读:

- `--qvalue` / `--matrix-qvalue`:蛋白/肽段鉴定与矩阵输出的 FDR 阈值
  (默认 1%);
- `--mass-acc`/`--mass-acc-ms1`:质量精度,须与仪器实际精度匹配,
  过松增假阳性、过严掉灵敏度;
- 固定/可变修饰决定搜谱空间(默认 Carbamidomethyl (C) 固定 +
  Oxidation (M) 可变,常规胶内/溶液酶切蛋白组);
- `diann_summary.py` 汇总每 run 鉴定数/定量数等 QC 指标
  (`diann_qc.tsv`)并生成蛋白 × 样本矩阵(`diann_protein_matrix.tsv`),
  供差异表达(如 MSstats/limma)与功能富集下游使用。

## 3. 参考文献

1. DIA-NN — https://github.com/vdemichev/DIA-NN
2. Demichev V, Messner CB, Volkmar B, et al. *DIA-NN: neural networks
   and interference correction enable deep proteome coverage in high
   throughput*. Nature Methods 17:1144-1150, 2020.
   DOI: 10.1038/s41592-020-00997-7.
