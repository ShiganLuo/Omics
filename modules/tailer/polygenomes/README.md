# tailer/polygenomes — 非 poly(A) RNA 3' 尾分析子模块

C 形态子模块(单规则 `tailer_global`):用 Tailer(jla-tailer)对 BAM 做
3' 端尾(tailing)分析,按 genome × sample 布局输出 tail 表。conda 环境
共享父目录 `../tailer.yaml`(见 `modules/tailer/`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}.bam` + `.bam.bai` | 比对 BAM/索引(形态 A 平铺布局,不随 genome 分层) |
| `genome.{genome}.gtf` | 每 genome 的 GTF 注释(`-a`;Tailer 据此推断成熟 3' 端,Ensembl GTF 效果良好) |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}_tail.csv   # Tailer tail 表
logdir/{genome}/{sample_id}/tailer.log               # 运行日志
```

Tailer 在输入 BAM 同目录产出 `{sample_id}_tail.csv`(运行前若有同名残留
文件会先删除),再 `shutil.copy2` 到声明输出;若工具未产出该文件则报
`FileNotFoundError`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.tailer`(`Tailer` 兼容键) | 可执行文件 | null(回退 `Tailer`) | Tailer 可执行路径 |
| `Params.tailer.read` | `-read` | 2 | 双端测序中携带 3' 端信息的是哪一端(1 或 2) |
| `Params.tailer.threshold` | `-t` | 100 | 距成熟 3' 端超过该距离(nt)的比对视为伪比对并丢弃 |
| `Params.tailer.rev_comp` | `--rev_comp` | false | 反向互补 reads(部分建库方法需要) |
| `genome.{genome}.gtf` | `-a` | — | GTF 注释(每 genome) |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录基路径 |

线程/内存经 `rule_res("tailer_global")` 配置:threads 默认 1,mem_mb 默认 8000。

引擎命令(每 genome × sample):

```bash
Tailer -a {genome.gtf} -read 2 -t 100 [ --rev_comp ] {indir}/{s}/{s}.bam
```

## 2. 算法原理与生物学意义

Tailer 面向非 poly(A) RNA 的 3' 端测序数据:在 global 模式下按 GTF 注释
推断各转录本的成熟 3' 端,把 read 3' 端与之比对,识别超出模板的 3' 尾
(如尿苷化、腺苷化等 post-transcriptional tailing),输出逐 read 的
tail 表(`{sample}_tail.csv`)。`-t` 阈值把偏离成熟端过远的比对判为伪
比对丢弃。生物学意义:非编码 RNA(snRNA、snoRNA、piRNA 等)的 3' 尾
修饰参与成熟、稳定与降解调控;tail 表可交给 Tailer-Analysis
(https://timnicholsonshaw.shinyapps.io/tailer-analysis/)做下游统计。

## 3. 参考文献

1. Tailer(3' end analysis of non-polyadenylated RNAs)—
   https://github.com/TimNicholsonShaw/tailer
2. jla-tailer(PyPI,`Tailer` 可执行来源)—
   https://pypi.org/project/jla-tailer/
3. Lykke-Andersen lab(UCSD)—
   https://labs.biology.ucsd.edu/lykkeandersen/index.html
