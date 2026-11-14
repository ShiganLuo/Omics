# scTE/polygenomes — 多基因组单细胞 TE 定量子模块(scTE)

polygenomes(每基因组独立目录)布局下的 scTE 封装:建 TE 索引
(`scTE_build`,资源型或下载型)→ 每样本每基因组从 Cell Ranger BAM
定量 TE 表达,输出 h5ad。共享父目录 `../scTE.yaml` 环境与
`modules/scTE/bin/scTE_quantify.py`;已知问题见父目录
`../KNOWN_ISSUES.md`(scTE 与 Cell Ranger 细胞数不一致等)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | Cell Ranger BAM 目录(上游 config 的 `outdir`):`{indir}/{sample_id}/{sample_id}.bam`(路径不含 genome 段) |
| `genome.<g>.gtf` / `.te_bed` (config) | 每基因组基因注释 GTF 与 TE 注释 BED(资源型建索引) |
| `genome.<g>.scTE_index` (config) | 可选预建索引;存在则直接用,否则用本模块构建的 `{outdir}/index/{Params.scTE.genome}.{mode}.idx` |
| `samples` / `outdir` / `logdir` / `logdir_index` (config) | 样本列表与目录 |

### 1.2 输出

```
{outdir}/index/{genome_name}.{mode}.idx                    # scTE 索引(命名按 Params.scTE.genome + mode)
{outdir}/{genome}/{sample_id}/{sample_id}_scTE.h5ad        # TE × 细胞表达 h5ad
```

回执:`{outdir}/index/scTE_build_index.sh`、
`{outdir}/{genome}/{sample_id}/scTE_quantify.sh`;日志
`{logdir_index}/scTE/scTE_build_index.log`、
`{logdir}/{genome}/{sample_id}/scTE_quantify.log`。

### 1.3 参数

config 位于 `Params.scTE`、`Procedure`、`genome` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.scTE.genome` | `-g`(下载型)/ 索引命名 | `hg38` | 基因组名(scTE_build 下载源) |
| `Params.scTE.mode` | `-m` | `exclusive` | 定量模式(scTE 语义:多映射读段处理策略,见 scTE 文档) |
| `Params.scTE.cb_tag` | `--cb-tag` | `CB` | 细胞 Barcode 标签 |
| `Params.scTE.umi_tag` | `--umi-tag` | `UB` | UMI 标签 |
| `genome.<g>.gtf` + `.te_bed` | `scTE_build -gene -te -g other` | — | 两者齐全走资源型建索引;否则 `scTE_build -g <genome>` 下载型 |
| `genome.<g>.scTE_index` | `--index` | 本模块构建的索引 | 预建索引优先 |
| `Procedure.scTE_build` | 索引引擎 | `scTE_build` | — |
| `Procedure.scTE` | `--scte-bin` | `scTE` | scTE 定量引擎(scTE_quantify.py 转调) |
| (threads) | `--threads` | 8(quantify)/ 1(index) | `rule_res` 可覆盖 |
| (mem_mb) | — | 16000(quantify)/ 8000(index) | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
# 资源型建索引(有 gtf + te_bed)
scTE_build -gene annotation.gtf -te TE.bed \
    -o out/index/hg38 -m exclusive -g other

# 下载型(仅有基因组名)
scTE_build -g hg38 -o out/index/hg38 -m exclusive

# 定量(父目录 bin/ 封装,内部转调 scTE)
python modules/scTE/bin/scTE_quantify.py \
    --input in/sample1/sample1.bam \
    --output out/GRCh38/sample1/sample1_scTE.h5ad \
    --index out/index/hg38.exclusive.idx \
    --cb-tag CB --umi-tag UB --threads 8 \
    --scte-bin scTE --sample-id sample1
```

注:`scTE_build_index` 的 params 通过 `wildcards.genome` 取
`config["genome"]` 的 gtf/te_bed,而该规则输出并无 genome wildcard
(索引按 `Params.scTE.genome` 命名)——接线时留意该引用的求值。

## 2. 算法原理与生物学意义

scTE 从单细胞 RNA-seq 比对中定量转座子(TE)家族表达:把读段按 TE
注释归属到 TE 家族/位点,以细胞 Barcode + UMI 计数,`exclusive` /
其他 mode 决定多重比对读段的分配策略(TE 高度重复,归属策略直接
影响定量保守/激进程度)。输出 h5ad(TE × 细胞)与 gene 矩阵互补,
用于单细胞水平的 TE 去抑制分析(如衰老、肿瘤、应激中 TE 家族异位
表达,以及 TE 转录本激活引发的先天免疫信号)。polygenomes 变体对
每个基因组各出一份 h5ad:不同物种/组装的 TE 注释差异大,结果须按
基因组分别解读。已知问题:scTE 判定的细胞数可能与 Cell Ranger
filtered 矩阵不一致(父目录 `KNOWN_ISSUES.md` 记录 GitHub issue #96
/#108 等),对细胞数敏感的下游需交叉核对。

## 3. 参考文献

1. scTE(单细胞 TE 定量工具)—— 模块文件内未附官方仓库链接;
   已知问题与其 GitHub issue 编号见 `modules/scTE/KNOWN_ISSUES.md`。
