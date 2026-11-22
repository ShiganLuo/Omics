# mimseq/tRNAtools — mim-tRNAseq 参考准备子模块(tRNA 解析/SNP 索引/比对索引)

mim-tRNAseq 流程的第一步:tRNA 序列与修饰位点解析、isodecoder 聚类、
修饰→SNP 索引生成、GSNAP 比对索引构建,状态序列化到 `state/` 供下游
(align/coverage/deseq)复用。引擎为共享父目录包装脚本
`modules/mimseq/bin/tRNAtools/run.py`(内部调用 vendored mimseq 包的
`tRNAtools` / `ssAlign`),conda 环境共享父目录 `../mimseq.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `Params.mimseq.species` | 内置物种名(如 `human`/`mouse`);非空时用 mimseq 自带参考数据,无需外部文件 |
| `genome.trnas` / `genome.trnaout`(config) | 自定义 tRNA 参考时的 tRNA 序列与注释文件(`species` 为空才消费) |
| `genome.mito_trnas` / `genome.plastid_trnas`(config) | 可选:线粒体/质体 tRNA 注释 |
| `Params.mimseq.data_dir` | mimseq 参考数据目录,缺省 `ROOT_DIR/modules/mimseq/mimseq/data` |
| — | 修饰位点数据 | 规则内固定 `ROOT_DIR/modules/mimseq/mimseq/modifications` |

### 1.2 输出

```
{outdir}/state/                      # tRNA 状态 pickle(聚类/索引等,下游消费)
{outdir}/mimseq_tRNAtools.sh         # 命令回执
```

收尾规则 `mimseq_tRNAtools_result` 消费 `state/`。

### 1.3 参数

config 位于 `Params.mimseq` 与 `genome` 下,与 `tRNAtools/run.py` CLI 对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `name` | `--name` | `tRNAseq` | 项目名 |
| `species` | `--species` | "" | 内置物种;空则用自定义 `--trnas/--trnaout` |
| `cluster_id` | `--cluster-id` | 0.97 | isodecoder 聚类相似度阈值 |
| `double_cca` | `--double-cca` | false | 双 CCA 末端处理开关 |
| `posttrans_mod_off` | `--posttrans` | false | 关闭转录后修饰建模 |
| `pretRNAs` | `--pretrnas` | false | 含 pre-tRNA |
| `local_modomics` | `--local-mod` | false | 使用本地 Modomics 修饰数据 |
| `snp_tolerance` | `--snp-tolerance` | true | 比对 SNP 容忍 |
| — | `--cluster` | true(`Params.mimseq.cluster`) | 聚类开关 |
| `Procedure.mimseq` | — | `mimseq` | mimseq 可执行路径(声明) |
| threads / mem_mb | — | 4 / 8000 | `rule_res("mimseq_tRNAtools")` 配置 |

注:解释器为规则 conda 环境(`../mimseq.yaml`)的 `python`,非 `Procedure` 配置。

直接调用示例(脱离流程):

```bash
python modules/mimseq/bin/tRNAtools/run.py \
  --modifications modules/mimseq/mimseq/modifications \
  --name tRNAseq --out out/ --species human \
  --data-dir modules/mimseq/mimseq/data \
  --threads 4 --cluster-id 0.97 --snp-tolerance --cluster
```

## 2. 算法原理与生物学意义

mim-tRNAseq 用一组特异引物捕获成熟 tRNA;由于 tRNA 修饰密集(逆转录
会在修饰位点停顿/错配),该流程把已知转录后修饰转化为参考序列上的
"修饰→SNP 索引",使修饰位点在比对中被容忍为系统错配,从而实现
isodecoder 分辨率的准确定量。tRNAtools 步骤即构建这套参考:解析 tRNA
序列/结构(`ssAlign`)、按 `cluster_id` 聚类 isodecoder、生成修饰 SNP
索引(`modsToSNPIndex`)与 GSNAP 索引(`generateGSNAPIndices`),并把
中间状态 pickle 到 `state/`。生物学意义:tRNA 表达与修饰动态关联翻译
调控、应激与疾病;isodecoder 级分辨率是传统 RNA-seq 做不到的。本模块
为流程封装。

## 3. 参考文献

1. mim-tRNAseq — https://github.com/nedialkova-lab/mim-tRNAseq
   (vendored 于 `modules/mimseq/mimseq/`,version.py 记 v1.3.6;
   作者 Drew Behrens, Nedialkova lab)。

## 已知问题 / 待办:`species` 未自动随基因组(暂缓修复)

`Params.mimseq.mimseq_tRNAtools.species` 目前是**静态配置值**,不随 meta/基因组自动变化
(对比:scanpy 的 `species` 已在 `scRNAseq.smk` 里传 `genome` 基因组键,自动跟随)。
本应按仓库范式(`GENOME_*_SPECIES` + `_fill_*`,如 `_fill_mei_rm_species`、VEP)从基因组
自动填充,但**暂缓**,原因(never-guess:不猜物种/取值):

1. **取值格式不确定**:mimseq 内置物种写作 `Hsap, Mmus, Scer`(mimseq/README.md),
   而本 README 上文示例写 `--species human`(human/mouse)——两者矛盾,需以 mimseq 实际
   `--species` 期望为准(查 argparse/文档)后才能填,填错会崩。
2. **物种支持范围有限**:mimseq 内置物种仅 Hsap/Mmus/Scer(人/鼠/酵母);恒河猴等其它物种
   需**自定义 tRNA 参考**(`species` 留空 + `--trnas`/`--trnaout`),此时不该自动填物种,
   而应走自定义参考路径。

**待办**:确认 mimseq `--species` 的准确取值(缩写 Hsap/Mmus 还是全名 human/mouse)后,
再按 `_fill_*` 范式接成"基因组→物种自动填充 + 缺映射则报错";非内置物种则保持自定义
tRNA 参考路径。在此之前请在配置里**显式**填 `species`(或留空用自定义 tRNA)。
