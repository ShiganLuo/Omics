# mimseq/polygenomes/tRNAtools — 多基因组 tRNA 参考与索引构建子模块(mimseq tRNAtools)

mim-tRNAseq 的 tRNAtools 步骤:解析 tRNA 序列、聚类同序列基因、生成 SNP
容忍索引与 GSNAP 比对索引,产出 `state/` 目录供 mimseq 主流程使用。
C 形态嵌套子模块,共享父目录 `../../mimseq.yaml` 环境、
`modules/mimseq/bin/tRNAtools/run.py` 脚本与 `modules/mimseq/bin/serialize.py`
(`from serialize import load_tRNA_state`,经 sys.path 引入)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `Params.mimseq.species` | 内置物种名(如 `mouse`);非空时无需外部 tRNA 文件(用 mimseq 内置数据) |
| `genome.<g>.trnas` / `genome.<g>.trnaout` | species 为空时的自定义 tRNA 序列与输出注释文件(必需) |
| `genome.<g>.mito_trnas` / `.plastid_trnas` | 线粒体/质体 tRNA 注释(可选) |

### 1.2 输出

```
{outdir}/{genome}/state/    # tRNAtools 状态目录(序列解析、SNP 索引、GSNAP 索引;directory())
```

日志:`{logdir}/{genome}/mimseq_tRNAtools.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录 |
| `Procedure.mimseq` | — | `mimseq` | mimseq 可执行路径(声明) |
| `Params.mimseq.species` | `--species` | `""` | 内置物种;空则用自定义 tRNA 文件 |
| `Params.mimseq.name` | `--name` | `tRNAseq` | 数据集名 |
| `Params.mimseq.data_dir` | `--data-dir` | `ROOT_DIR/modules/mimseq/mimseq/data` | mimseq 数据目录 |
| `Params.mimseq.cluster_id` | `--cluster-id` | `0.97` | 序列聚类一致性阈值 |
| `Params.mimseq.snp_tolerance` | `--snp-tolerance` | `true` | 默认传 |
| `Params.mimseq.cluster` | `--cluster` | `true` | 默认传 |
| `Params.mimseq.double_cca` | `--double-cca` | `false` | true 时传 |
| `Params.mimseq.posttrans_mod_off` | `--posttrans` | `false` | true 时传 |
| `Params.mimseq.pretRNAs` | `--pretrnas` | `false` | true 时传 |
| `Params.mimseq.local_modomics` | `--local-mod` | `false` | true 时传 |
| `genome.<g>.trnas` / `trnaout` | `--trnas` / `--trnaout` | — | 自定义模式输入 |
| `genome.<g>.mito_trnas` / `plastid_trnas` | `--mito-trnas` / `--plastid-trnas` | `""` | 非空才传 |
| 固定 | `--modifications ROOT_DIR/modules/mimseq/mimseq/modifications`、`--threads` | — | 修饰知识库目录 |
| 线程/内存 | `--threads` | `threads=4`,`mem_mb=8000` | `rule_res(...)` 可覆盖 |

环境:共享父目录 yaml `../../mimseq.yaml`;run 块内固定 `source
/home/zhousha/miniforge3/... /data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_`
后以其 python 运行(环境现状,迁移需替换)。

直接调用示例(等价引擎命令):

```bash
python modules/mimseq/bin/tRNAtools/run.py \
  --modifications modules/mimseq/mimseq/modifications \
  --name tRNAseq --out output/GRCh38/ --species mouse \
  --data-dir modules/mimseq/mimseq/data --threads 4 --cluster-id 0.97 \
  --snp-tolerance --cluster [--double-cca] [--posttrans] [--pretrnas] [--local-mod]
```

## 2. 算法原理与生物学意义

纯工具封装(mim-tRNAseq 的 tRNAtools):tRNA 基因家族高度同源且多拷贝,
tRNAtools 先解析 tRNA 序列(含 intron、CCA 尾),按 `cluster_id` 聚类同序列
基因形成唯一比对单元,再建容忍 SNP/修饰位点的比对索引,使测序读段能正确
归属到 tRNA 基因/序列变体。`state/` 目录是下游定量与修饰分析的共享参考。
结果怎么读:`state/` 内为序列与索引文件,无需逐文件解读;确认目录生成即
参考就绪。

## 3. 参考文献

1. mimseq / mim-tRNAseq — https://github.com/goodarzilab/mimseq
2. GMAP/GSNAP — Wu TD, Nacu S. *Fast and SNP-tolerant detection of complex
   variants and splicing in short reads*. Bioinformatics 26(7):873–881, 2010.
   DOI: 10.1093/bioinformatics/btq057
