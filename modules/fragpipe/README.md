# fragpipe（MSFragger + Philosopher：FragPipe 生态搜库与 FDR）

按 PXD053575 / Blood 2025;147(23):2809 原文方法重跑鉴定：MSFragger 搜库
（质量校准、TMTpro 修饰、strict trypsin、漏切 2）→ Philosopher
PeptideProphet + picked protein FDR 过滤 → psm/peptide/protein 报表。

两个子工具共用 `fragpipe.yaml` 环境（SIF: `/data/16T/Database/env/fragpipe/fragpipe.sif`）。

## 规则

| 规则 | 输出 |
|---|---|
| msfragger_run | `{outdir}/{sample}/{sample}.pepXML` + `{sample}.tsv`（每样本） |
| msfragger_result | 聚合目标 |
| philosopher_database | `{outdir}/database/fragpipe_db.fasta`（target + rev_ decoy + contaminants） |
| philosopher_filter | `{outdir}/{sample}/interact-{sample}.pep.xml` + `psm.tsv` + `peptide.tsv`（+ `protein.tsv`/`ion.tsv`） |
| philosopher_result | 聚合目标 |

## 方法与参数（来自原始 pepXML 的参数转储）

`bin/make_fragger_params.py` 生成 MSFragger 参数文件，默认值即原文值：

- 质量：precursor -20/+20 ppm、CORRECTED（`calibrate_mass=1`）、
  `precursor_true_tolerance 10`、isotope_error -1/0/1/2、deisotope
- 碎片：0.02 Da、b,y 离子、max_fragment_charge 1
- 酶：stricttrypsin（KR|P）、`num_enzyme_termini 2`、漏切 2
- 固定修饰：Carbamidomethyl(C) 57.021464、**TMTpro(K) 304.207**
- 可变修饰：Ox(M) 3、protein N-term 乙酰化 42.010565、
  **TMTpro 任意 N 端 304.207**、pyro-Glu(-17 nQnC / -18 nE)
- 其它：digest_mass_range 350-10000、precursor_charge 2-6、
  delta_mass_exclude_ranges (-1.5,3.5)、max_variable_mods 5
- 输出：`tsv_pepxml_pin`

原文用 MSFragger 3.2/3.5；bioconda 提供 4.4.1（参数兼容，已注明版本差异）。

## config

```
genome.search_fasta        target 蛋白 FASTA（不含 decoy；philosopher 生成 rev_ decoy + contam）
Params.msfragger.threads   默认 16
Params.msfragger.label     tmt16（默认，原文方法）/ label-free（去 TMT 修饰）
Params.msfragger.jar       MSFragger-3.x jar 路径（推荐：与原文版本一致且无需
                           license key，经 java -jar 直接运行）
Params.msfragger.license_key  仅 4.x wrapper 模式必填（免费学术授权：
                           msfragger-upgrader.nesvilab.org/upgrader/）
Params.msfragger.java_xmx  jar 模式 JVM 内存（默认 16G）
Params.philosopher.*       pep_fdr / prot_fdr / ion_fdr（0.01）/ picked（true）
```

## 填充清单

| 键 | 文件 | 来源 |
|---|---|---|
| `genome.search_fasta` | 蛋白库 | 原文为 GENCODE hsV42 basic + contam + entrezgene decoys；本地可用 `/data/16T/Database/Reference/human/GENCODE/GRCh38/gencode.v49.pc_translations.fa`（philosopher 自动补 contam+decoy）。若能拿到核心的原始 .fas 更佳 |
| `Params.msfragger.jar` 或 `license_key` | jar / 字符串 | **推荐**：从蛋白组核心拷贝 `MSFragger-3.5.jar`（他们就是用它搜的这批数据；jar 仅需 Java，3.x 无 key 要求）。否则在 msfragger-upgrader.nesvilab.org/upgrader/ 免费注册拿 key 用 4.4.1 |

## 运行示例（独立跑）

```bash
snakemake -s /data/16T/luosg/workflow/Omics/modules/fragpipe/msfragger/msfragger.smk \
    --configfile msfragger_run.json --cores 16 --sdm apptainer \
    --singularity-args "--bind /data/16T" --rerun-triggers=mtime msfragger_result

snakemake -s /data/16T/luosg/workflow/Omics/modules/fragpipe/philosopher/philosopher.smk \
    --configfile philosopher_run.json --cores 4 --sdm apptainer \
    --singularity-args "--bind /data/16T" --rerun-triggers=mtime philosopher_result
```

独立运行的 configfile 需携带 `ROOT_DIR` 与 `env`（见 modules.md）。

## 注意

- MSFragger 把结果写在输入文件旁：规则用硬链接把 mzML 引入样本工作目录，
  输出绝不污染 2_raw2mzml。
- Philosopher 的 `database` 输出文件名带日期，规则生成后即改名固定；
  `database --annotate` 必须在每个样本的 workspace 里执行（filter 的蛋白层
  读取本地 `.meta/db.bin`）。
- 玩具 FASTA 头会触发 Philosopher 5.1.2 的切片越界 panic
  （UpdateLayerswithDatabase）——蛋白头须是 GENCODE/UniProt 格式。
