# telomere — 端粒长度分析模块(Telogator2 / 组装扫描 / 读段密度 / tidk)

四种端粒长度测量方法并存:① Telogator2 逐染色体臂端粒长度;② 组装 contig
末端端粒重复扫描(推荐小鼠 30–150 kb 端粒,HiFi 读段 15–25 kb 跨不全);
③ 读段级 k-mer 密度全基因组平均;④ tidk 社区工具组装扫描。模块私有脚本
`bin/scan_assembly_telomere.py`、`bin/read_density_telomere.py`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_substring` | 形态 A 拼接:`{indir}/{sample_id}/{sample_id}[.{bam_substring}].bam`(+ `.bai`,Telogator2 需要) |
| `assembly_dir` | hifiasm 组装目录(centromere 模块产出);缺省回退 `{outdir}/../centromere/{sample_id}/assembly/asm.bp.p_ctg.fa` |
| `samples` | 样本列表 |
| tidk 数据库 | `~/.local/share/tidk/tidk_database.csv`;不存在时先跑 `tidk_init`(`tidk build`)建库 |

### 1.2 输出

```
{outdir}/{sample_id}/telogator2/tlens_by_allele.tsv             # 逐等位端粒长度
{outdir}/{sample_id}/telogator2/all_final_alleles.png           # 等位图
{outdir}/{sample_id}/telogator2/violin_atl.png                  # 端粒长度小提琴图
{outdir}/{sample_id}/assembly_scan/{sample_id}_assembly_telomere.tsv        # contig 末端端粒
{outdir}/{sample_id}/assembly_scan/{sample_id}_assembly_telomere_stats.txt  # 汇总
{outdir}/{sample_id}/read_density/{sample_id}_read_telomere.tsv
{outdir}/{sample_id}/read_density/{sample_id}_read_telomere_stats.txt
{outdir}/{sample_id}/tidk/{sample_id}_tidk_telomeres.tsv        # tidk 结果(原生 tsv 改名)
{outdir}/.tidk_build.done                                       # tidk 建库标记
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `logs` | 目录 |
| `bam_substring` | 文件名插入段 | `""` | 见输入 |
| `assembly_dir` | — | `""` | 见输入 |
| `Params.telogator2.species` | `-t` 参考选择 | `human` | 非 human 时自动找 telogator2 安装内 `non-human/telogator-ref-<species>.fa.gz`,找不到回退默认 |
| `Params.telogator2.assembly_scan_length` | `--scan_length` | `50000` | 组装末端扫描长度(bp) |
| `Params.telogator2.n_chrom_arms` | `--n_chrom_arms` | `40` | 染色体臂数(读段密度法分母) |
| `Params.telogator2.threads` | — | 16(rule_res) | Telogator2 线程 |
| `Procedure.telogator2` / `tidk` | `telogator2 -i` / `tidk search` | — | 可执行路径 |
| 固定 | `tidk search --string TTAGGG` | — | 端粒 motif |

脚本 CLI 补充:`scan_assembly_telomere.py` 另有 `--window_size`、`--step_size`、
`--min_tel_frac`;`read_density_telomere.py` 另有 `--min_read_length`、
`--min_tel_kmers`。环境:`conda: "telomere.yaml"`(模块目录内)。

直接调用示例(等价引擎命令 + 模块脚本):

```bash
telogator2 -i s1.bam                                  # 非人:加 -t telogator-ref-<species>.fa.gz
python modules/telomere/bin/scan_assembly_telomere.py \
  --fasta asm.bp.p_ctg.fa --sample_name s1 --output_dir out/s1/assembly_scan \
  --scan_length 50000
python modules/telomere/bin/read_density_telomere.py \
  --bam s1.bam --sample_name s1 --output_dir out/s1/read_density \
  --n_chrom_arms 40
tidk build && tidk search --string TTAGGG --output s1_tidk_search \
  --dir out/s1/tidk asm.bp.p_ctg.fa
```

## 2. 算法原理与生物学意义

端粒(TTAGGG 重复)长度是细胞衰老与基因组稳定性的经典指标;小鼠端粒
长达 30–150 kb,超出 HiFi 读长,故模块并行四种方法:

- **Telogator2**:读段种子-延伸匹配端粒/亚端粒参考,逐染色体臂、逐等位
  估计长度(需读段跨过端粒-亚端粒交界,适合人端粒);
- **Approach A 组装扫描**:在 contig 两端 `scan_length` 窗口内找连续
  CCCTAA/TTAGGG 重复段,给每 contig 每端端粒长(小鼠首选);
- **Approach B 读段密度**:全读段 telomeric k-mer 总碱基数 ÷ 染色体臂数
  → 全基因组平均端粒长;
- **Approach C tidk**:社区标准工具在组装上搜端粒重复(`tidk search`)。

结果怎么读:`*_stats.txt` 给平均端粒长;`tlens_by_allele.tsv` 逐等位长度
(端粒等位间差异大);A/B/C 互为验证,小鼠以 A 为准。

## 3. 参考文献

1. Telogator2 — https://github.com/telomere-toolkit(Telomere toolkit;
   如引用请以官方页面为准)
2. tidk — https://github.com/tolkit/tidk
3. 端粒与衰老综述:López-Otín C et al. *Hallmarks of aging: An expanding
   universe*. Cell 186(2):243–278, 2023. DOI: 10.1016/j.cell.2022.11.001
