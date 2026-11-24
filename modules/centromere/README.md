# centromere — 着丝粒卫星 DNA 分析模块(RepeatMasker + 统计)

对上游 hifiasm 组装(`modules/hifiasm`)的 contig 用 RepeatMasker 注释卫星
重复,再提取着丝粒卫星统计。两步:RepeatMasker 重复注释 →
`bin/extract_centromere_stats.py` 统计。组装不属本模块:下游接线经
`assembly_dir` 配置键(或共享 outdir)指向 `modules/hifiasm` 的产出布局。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `assembly_dir` | 组装目录(可选,缺省 = `outdir`):期望布局 `{assembly_dir}/{sample_id}/assembly/asm.bp.p_ctg.fa`,即 `modules/hifiasm` 的产出契约;subworkflow 用模块级常量引用上游 config 的 `outdir` |
| `samples` | 样本列表 |
| `Params.RepeatMasker.dfam_h5` | 可选:Dfam FamDB(h5)库;给了则 `repeatmasker_init` 先建库缓存(链接 h5 进 RepeatMasker famdb 库目录、预热缓存,避免并发 makeblastdb 冲突),`repeatmasker_run` 依赖该标记;不给则直接用 RepeatMasker 自带默认库 |

### 1.2 输出

```
{outdir}/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.out              # 重复注释明细
{outdir}/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.tbl              # RepeatMasker 汇总表
{outdir}/{sample_id}/{sample_id}.centromere_stats.txt              # 着丝粒卫星统计
{outdir}/.repeatmasker_lib_init.done                               # 库缓存完成标记(配置 dfam_h5 时)
```

日志:单样本规则 `{logdir}/{sample_id}/{规则名}.log`
(`repeatmasker_run.log` / `centromere_extract.log`);
多样本规则 `repeatmasker_init` 走 `logdir_combine + "/centromere/repeatmasker_init.log"`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `outdir` / `logdir` / `logdir_combine` | — | `output` / `logs` / 同 `logdir` | 目录;混合模块两个日志根都要传 |
| `assembly_dir` | — | 同 `outdir` | 组装 FASTA 目录,见输入 |
| `Params.RepeatMasker.species` | `-species` / `--species` | 无(必填,禁止静默猜测) | 重复库物种(如 `Mus musculus`);同时驱动 `extract_centromere_stats.py` 的着丝粒家族表(接受 `mouse`/`human` 别名;其他物种降级为通用卫星统计)。node.py 可按 genome 自动填充 |
| `Params.RepeatMasker.repeatmasker_init.species` / `.repeatmasker_run.species` | 同上 | 回退到 `.species` | 按规则覆盖物种(init 与 run 可各自指定) |
| `Params.RepeatMasker.dfam_h5` | `--dfam-h5` | `null` | Dfam FamDB 库,见输入 |
| `Procedure.RepeatMasker` / `python` | — | 同名命令 | 可执行路径 |
| 线程/内存 | — | RepeatMasker 24/`mem_mb=32000` | 真实入口是 `resources.<规则名>.threads/.mem_mb`(`rule_res()` 解析并按宿主预算钳制) |

环境:`conda: "centromere.yaml"`(模块目录内)。

直接调用示例(模块引擎命令,与回执脚本一致):

```bash
# RepeatMasker 库缓存初始化(配置 dfam_h5 时先跑一次)
python modules/centromere/bin/repeatmasker_init.py \
  --repeatmasker RepeatMasker --species 'Mus musculus' --threads 1 \
  --dfam-h5 /path/dfam39_full.7.h5 --marker out/.repeatmasker_lib_init.done

# RepeatMasker 重复注释(单命令,直接调用工具)
RepeatMasker -species 'Mus musculus' -pa 24 -dir out/s1/RepeatMasker \
  -gff out/s1/assembly/asm.bp.p_ctg.fa

# 着丝粒卫星统计(--species 接受 mouse/human 或 'Mus musculus'/'Homo sapiens';
# 无专属家族表的物种(如 'Macaca mulatta')降级为通用卫星统计。
# 另有 --max_gap 50000 --min_block_len 100000)
python modules/centromere/bin/extract_centromere_stats.py \
  --rm_out out/s1/RepeatMasker/asm.bp.p_ctg.fa.out \
  --output out/s1/s1.centromere_stats.txt --species 'Mus musculus'
```

## 2. 算法原理与生物学意义

着丝粒由高阶卫星重复(小鼠 Major/Minor satellite,人类 α-satellite)构成,
短读长无法解析,需 HiFi 长读长组装(上游 `modules/hifiasm`:hifiasm 串图
(BOG)组装 HiFi 读段,取主 contig `bp.p_ctg`)。RepeatMasker 按物种重复库
(如 Dfam)注释卫星 DNA;`extract_centromere_stats.py` 按物种着丝粒家族表
(如小鼠 MajorSatellite: MSAT/Major_satellite/GSAT…)把散落注释聚成块
(`--max_gap` 缺口合并、`--min_block_len` 最小块长),输出着丝粒卫星的
长度/组成统计。结果怎么读:`centromere_stats.txt` 为每 contig 的卫星块
统计,可看各样本着丝粒大小/卫星组成差异(着丝粒重复拷贝数与染色体稳定
性、衰老相关异染色质变化相关)。

## 3. 参考文献

1. RepeatMasker — https://www.repeatmasker.org/ ,
   https://github.com/rmhubley/RepeatMasker
2. Dfam — https://www.dfam.org/(Storer J et al.,
   Nucleic Acids Research 49(D1):D80–D87, 2021)
3. hifiasm — https://github.com/chhylp123/hifiasm(Cheng H et al.,
   Nature Methods 18:608–610, 2021. DOI: 10.1038/s41592-020-01056-5)
