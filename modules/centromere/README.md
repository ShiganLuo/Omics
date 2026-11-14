# centromere — 着丝粒卫星 DNA 分析模块(hifiasm + RepeatMasker)

对 HiFi BAM 从头组装,用 RepeatMasker 注释卫星重复,再提取着丝粒卫星
统计。三步:hifiasm 组装 → RepeatMasker 重复注释 →
`bin/extract_centromere_stats.py` 统计。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_substring` | 形态 A 拼接:`{indir}/{sample_id}[.{bam_substring}].bam` + 同名 `.bam.pbi` 索引(PacBio HiFi BAM) |
| `samples` | 样本列表 |
| `Params.RepeatMasker.dfam_h5` | 可选:Dfam FamDB(h5)库;给了则先跑 `repeatmasker_init` 建库缓存(避免并发 makeblastdb 冲突) |

### 1.2 输出

```
{outdir}/{sample_id}/assembly/asm.bp.p_ctg.fa                      # 组装主 contig
{outdir}/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.out              # 重复注释明细
{outdir}/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.tbl              # RepeatMasker 汇总表
{outdir}/{sample_id}/{sample_id}.centromere_stats.txt              # 着丝粒卫星统计
{outdir}/.repeatmasker_lib_init.done                               # 库缓存完成标记
```

日志:`{logdir}/{sample_id}/hifiasm.log`、`repeatmasker.log`、
`centromere_extract.log`、`{logdir}/all/RepeatMasker/repeatmasker_init.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `logs` | 目录 |
| `bam_substring` | 文件名插入段 | `""` | 见输入 |
| `Params.hifiasm.threads` | `-t` | 12(rule_res) | hifiasm 线程 |
| `Params.RepeatMasker.species` | `-species` | `Mus musculus`(json 示例 `mouse`) | 重复库物种 |
| `Params.RepeatMasker.threads` | `-pa` | 12(rule_res) | RepeatMasker 并行 |
| `Params.RepeatMasker.dfam_h5` | — | `""` | Dfam FamDB 库,见输入 |
| `Procedure.hifiasm` / `RepeatMasker` / `python` | — | `hifiasm` / `RepeatMasker` / `python` | 可执行路径 |
| 线程/内存 | — | hifiasm `mem_mb=64000`;RepeatMasker `mem_mb=32000` | `rule_res(...)` 可覆盖 |

环境:`conda: "centromere.yaml"`(模块目录内)。

直接调用示例(等价引擎命令 + 模块脚本):

```bash
samtools fastq -@ 12 s1.bam | gzip -c > s1.fq.gz
hifiasm -o asm -t 12 s1.fq.gz
awk '/^S/{print ">"$2; print $3}' asm.bp.p_ctg.gfa > asm.bp.p_ctg.fa
RepeatMasker -species 'Mus musculus' -pa 12 -dir out/s1/RepeatMasker -gff asm.bp.p_ctg.fa
python modules/centromere/bin/extract_centromere_stats.py \
  --rm_out out/s1/RepeatMasker/asm.bp.p_ctg.fa.out \
  --output out/s1/s1.centromere_stats.txt   # 另有 --species mouse --max_gap 50000 --min_block_len 100000
```

## 2. 算法原理与生物学意义

着丝粒由高阶卫星重复(小鼠 Major/Minor satellite,人类 α-satellite)构成,
短读长无法解析,需 HiFi 长读长组装。hifiasm 用串图(BOG)组装 HiFi 读段,
取主 contig(`bp.p_ctg`);RepeatMasker 按物种重复库(如 Dfam)注释卫星
DNA;`extract_centromere_stats.py` 按物种着丝粒家族表(如小鼠
MajorSatellite: MSAT/Major_satellite/GSAT…)把散落注释聚成块
(`--max_gap` 缺口合并、`--min_block_len` 最小块长),输出着丝粒卫星的
长度/组成统计。结果怎么读:`centromere_stats.txt` 为每 contig 的卫星块
统计,可看各样本着丝粒大小/卫星组成差异(着丝粒重复拷贝数与染色体稳定
性、衰老相关异染色质变化相关)。

## 3. 参考文献

1. hifiasm — https://github.com/chhylp123/hifiasm(Cheng H et al.,
   Nature Methods 18:608–610, 2021. DOI: 10.1038/s41592-020-01056-5)
2. RepeatMasker — https://www.repeatmasker.org/ ,
   https://github.com/rmhubley/RepeatMasker
3. Dfam — https://www.dfam.org/(Storer J et al.,
   Nucleic Acids Research 49(D1):D80–D87, 2021)
