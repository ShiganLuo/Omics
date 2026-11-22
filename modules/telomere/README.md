# telomere — 端粒分析父模块(family)

端粒分析家族,按 `modules/modules.md` layout C 组织:父目录 `telomere/` 共享 env(`telomere.yaml`)
与 helper,下设子模块。输出统一在 **`{outdir}/telomere/{sub}/{method}/{sample}/`**(method 在样本前)。

## 子模块

| 子模块 | 作用 | 输出 |
|---|---|---|
| [`length/`](length/) | **端粒长度**:telogator2 / assembly_scan / read_density / tidk + 深度 QC + 合并图 | `{outdir}/telomere/length/{method}/{sample}/…` |
| [`composition/`](composition/) | **端粒组成/变异基序(TVR/ITS)+ 长度**(telomere_profile,多样本) | `{outdir}/telomere/composition/…` |
| [`regions/`](regions/) | **端粒区域 BED**(telomere_core/subtelomere,参考扫描) | `{outdir}/telomere/regions/{sample}/…` |

## 共享约定

- **env/自含**:`conda: "../telomere.yaml"`、`container: sif("../telomere.yaml", …)`;
  子模块 `.smk` `include: "../../common/common.smk"`;`bin/` 自包含(规范#0)。
- **输入(字典模式)**:`inputs[<sub>][<key>][<sample>]`(仿 sexcheck `sample_bams`),
  缺失回退 shape A(`indir`+`bam_substring`)。见 `subworkflow/subworkflow.md`「传参规范」。
- **参数**:`Params.telomere` 公共字段(`enabled`/`method`/`formats`)在顶,方法/子模块字段
  (如 `telogator2.species`、`read_density.n_chrom_arms`)在各自键下。
- **method**:`Params.telomere.method` 默认**全部方法**;缺输入的方法按 PacVar 接线约定处理。
- **绘图**:每方法一张**全样本合并图**(`combined.{fmt}`);4 法全指定时加**跨方法总览**
  (`combined_methods.{fmt}`)。多格式 `Params.telomere.formats`(规范#9)。
- **深度 QC(⑥)**:输入比对 BAM 的方法用**参考端粒坐标 + 测序深度**作质控(见 `length/`)。

## 方法选择

模块定义全部 rule、不写条件分支;`Params.telomere.method`(列表)由 `node.py` 注册输出决定
跑哪些(「输出决定流程走向」)。详细方法原理、结果解读、图像元件见各子模块 `README.md`。
