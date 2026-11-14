# subsample — 高丰度小 RNA 样本降采样模块(seqtk)

按样本名识别高丰度小 RNA 样本,用 `seqtk sample` 把其双端 reads 降采样到
设定条数;其余样本符号链接直通(单规则 `subsample_fastq`)。

## 1. 模块输入输出与参数

### 1.1 输入

形态 A:`indir` + 样本列表,仅双端:

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}_1.fq.gz` | R1,gzip 压缩 FASTQ |
| `{indir}/{sample_id}/{sample_id}_2.fq.gz` | R2 |

判定规则:样本 ID(大写)包含任一 `abund_small_rnas` 名称(大写)即视为
高丰度样本(`is_abundant` 子串匹配)。

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_1.fq.gz   # 降采样 R1,或指向输入的符号链接
{outdir}/{sample_id}/{sample_id}_2.fq.gz   # 降采样 R2,或符号链接
logdir/{sample_id}/subsample.log           # 运行日志
```

高丰度样本走 zcat → `seqtk sample` → gzip;非高丰度样本直接
`os.symlink` 到输入文件,零拷贝通过。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.subsample.abund_small_rnas` | — | [] | 高丰度小 RNA 名称列表(样本名子串匹配,大小写不敏感) |
| `Params.subsample.subsample_n` | `seqtk sample` 的条数参数 | 100000 | 高丰度样本保留 reads 数;非高丰度样本传 0(全量直通) |
| `Params.subsample.subsample_seed` | `-s` | 42 | 随机种子;R1/R2 同种子同条数以保持配对 |
| `indir` | — | `input` | 输入目录 |
| `outdir` | — | `output` | 输出目录 |
| `logdir` | — | `logs` | 日志目录 |

注意:`subsample.json` 中另含 `abundant_rnas` / `n` / `seed` / `direction` /
`cut_length` 键,**与 .smk 读取的键名不一致且 .smk 未读取**;实际生效键以
上表为准。

线程/内存经 `rule_res("subsample_fastq")` 配置:threads 默认 1,mem_mb 默认 8000。

引擎命令(高丰度样本):

```bash
seqtk sample -s 42 {s}_1.tmp.fq 100000 > {s}_1.sub.fq
seqtk sample -s 42 {s}_2.tmp.fq 100000 > {s}_2.sub.fq
gzip -c {s}_1.sub.fq > {outdir}/{s}/{s}_1.fq.gz
gzip -c {s}_2.sub.fq > {outdir}/{s}/{s}_2.fq.gz
```

## 2. 算法原理与生物学意义

纯工具封装:`seqtk sample` 用固定种子做随机抽取,同种子作用于 R1/R2 时
选中相同的 read 编号集合,双端配对关系保持完整。生物学意义:小 RNA 测序
中个别物种(如 rRNA、tRNA 片段)reads 可占绝大部分,降采样高丰度样本能
压缩测序量/计算量、平衡样本间深度,而低丰度样本保留全量避免丢信息。

## 3. 参考文献

1. seqtk — https://github.com/lh3/seqtk
