# decoil — ecDNA 重建与可视化模块(DeCoil + 自研绘图)

从 PacBio HiFi 长读长重建染色体外环状 DNA(ecDNA/eccDNA):DeCoil
`reconstruct` 消费 BAM + 外部 SV VCF + 覆盖 bigwig 产出重建环,自研
`bin/` 脚本做统计图、跨样本比较、环结构图与 IGV/UCSC 浏览。**本模块
只做重建**:SV calling(pbsv/sniffles)与覆盖度轨道须由上游提供
(推荐 pbsv → VCF、igv.wig/deeptools → bigwig)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_substring` (config) | 坐标排序 BAM + 索引:`{indir}/{sid}/{sid}{.bam_substring}.bam` 与同名 `.bai`(形态 A;示例 `bam_substring=sorted_markdup`) |
| `vcf_dir` + `vcf_suffix` (config) | 外部 SV VCF:`{vcf_dir}/{sid}/{sid}{vcf_suffix}`(默认 `.sv.vcf.gz`;pbsv 产物先经 `bin/pbsv2decoil.py` 转成 sniffles1 风格 plain VCF——DeCoil 需要 STRANDS/CHR2/END 与 FORMAT DR/DV 字段,且不读 bgzip) |
| `bw_dir` + `bw_suffix` (config) | 覆盖 bigwig:`{bw_dir}/{sid}/{sid}{bw_suffix}`(`bw_dir` 缺省 `outdir`,默认 `.bigwig`) |
| `genome.fasta` / `genome.gtf` (config) | 参考 FASTA 与基因注释;**两者必填且须存在**(load 期 `ValueError`) |
| `Params.decoil.igv` (config) | IGV 配置段(fastaURL/indexURL/GTF 轨道 + publicPathMap,参照 modules/track/track.json);`decoil_igv` 缺它直接 `ValueError` |
| `samples` / `outdir` / `logdir` / `ROOT_DIR` (config) | 样本列表与目录 |

### 1.2 输出

```
# decoil_reconstruct(每样本)
{outdir}/{sid}/decoil/reconstruct.bed                    # 全部重建片段(按序)
{outdir}/{sid}/decoil/reconstruct.ecDNA.bed              # ecDNA 标记环的片段
{outdir}/{sid}/decoil/reconstruct.ecDNA.filtered.bed     # 过 --filter-score 的环
{outdir}/{sid}/decoil/summary.txt                        # 每环一行:topology + proportion
{outdir}/{sid}/decoil/decoil_input.vcf                   # pbsv2decoil 转换后的输入 VCF

# decoil_plot / decoil_structure / decoil_igv / decoil_ucsc(每样本)
{outdir}/{sid}/visualization/decoil_stats_ecdna_regions.tsv + decoil_stats.<fmt>
{outdir}/{sid}/visualization/decoil_structure_genes.tsv + decoil_structure.<fmt>(可每环一文件)
{outdir}/{sid}/visualization/igv_ecdna.html(+ igv_ecdna_config.json)
{outdir}/{sid}/visualization/ucsc_track.txt

# decoil_compare(全样本)
{outdir}/compare/decoil_compare.png
{outdir}/compare/decoil_compare_shared_unique.tsv        # 共有/独有区段
{outdir}/compare/decoil_compare_ecdna_regions.tsv
```

`<fmt>` 按各段 `formats` 展开(默认 png)。回执:各规则 `*.sh` 落在对应
输出目录;日志 `{logdir}/{sid}/decoil_*.log`、`{logdir}/decoil_compare.log`。

### 1.3 参数

config 位于 `Params.decoil` 下,规则级键(reconstruct/plot/compare/
structure/igv)优先,回落到父级同名键(如 `formats`):

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.decoil.reconstruct.sv_caller` | `--sv-caller` | `sniffles1` | 输入 VCF 的来源格式(sniffles1/sniffles2/lumpy/multivcf) |
| `...min_sv_len` | `--min-sv-len` | 500 | 最小 SV 长度(bp) |
| `...fragment_min_cov` / `fragment_max_cov` | `--fragment-min-cov` / `--fragment-max-cov` | 5 / 100000 | 片段覆盖度上下限 |
| `...fragment_min_size` | `--fragment-min-size` | 500 | 最小片段长度(bp) |
| `...min_vaf` | `--min-vaf` | 0.01 | 最小等位频率 |
| `...min_cov_alt` / `...min_cov` | `--min-cov-alt` / `--min-cov` | 6 / 8 | alt / 总覆盖度下限 |
| `...max_explog_threshold` | `--max-explog-threshold` | 0.1 | exp-log 阈值 |
| `...filter_score` | `--filter-score` | 5 | `filtered.bed` 的环得分阈值 |
| `...extend_allowed_chromosomes` | `--extend-allowed-chr` | "" | 非空才传 |
| `...skip` | `--skip` | false | 跳过 FASTA 相关步骤 |
| `Params.decoil.plot.formats` | `decoil_plot.py -f` | `[png]` | 统计图格式(可多) |
| `Params.decoil.compare.formats` | `decoil_plot.py -f` | `[png]` | 比较图格式 |
| `Params.decoil.structure.formats` | `decoil_structure.py -f` | `[png]` | 结构图格式;**decoil_plot 的 output 文件列表也取本键**(若与 plot.formats 不一致会与实际产物不符,建议保持一致) |
| `Params.decoil.structure.max_per_fig` | `--max-per-fig` | 12 | 每图最多环数 |
| `Params.decoil.structure.one_file_per_cycle` | `--one-file-per-cycle` | true | 每环单独出图 |
| `Params.decoil.igv.locus` | `--locus` | "" | IGV 初始定位;缺省打开最大 ecDNA 区段 |
| `Procedure.decoil` | 引擎可执行 | `decoil` | DeCoil 2.x(无 --threads 参数) |
| (threads) | — | 16(reconstruct)/ 4(plot/structure)/ 1/2 | `rule_res` 可覆盖 |
| (mem_mb) | — | 32000(reconstruct)/ 12000/8000/4000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
# pbsv VCF → sniffles1 风格(decoil reconstruct 不能直接读 pbsv)
python3 modules/decoil/bin/pbsv2decoil.py sample.sv.vcf.gz decoil_input.vcf

decoil reconstruct -b sample.sorted_markdup.bam -i decoil_input.vcf \
    -c sample.bigwig -r genome.fa -g annotation.gtf \
    -o out/sample/decoil --name sample \
    --sv-caller sniffles1 --min-sv-len 500 \
    --fragment-min-cov 5 --fragment-max-cov 100000 --fragment-min-size 500 \
    --min-vaf 0.01 --min-cov-alt 6 --min-cov 8 \
    --max-explog-threshold 0.1 --filter-score 5

# 统计图 / 比较图(同一脚本两种 mode)
python modules/decoil/bin/decoil_plot.py --mode single \
    --sample s1:out/s1/decoil/summary.txt \
    --ecdna-bed s1:out/s1/decoil/reconstruct.ecDNA.filtered.bed \
    -o out/s1/visualization/decoil_stats -f png

# 环结构图
python modules/decoil/bin/decoil_structure.py --sample s1 \
    --ecdna-bed out/s1/decoil/reconstruct.ecDNA.filtered.bed \
    --summary out/s1/decoil/summary.txt --gtf annotation.gtf \
    --bw out/s1/s1.bigwig \
    -o out/s1/visualization/decoil_structure --max-per-fig 12 \
    --one-file-per-cycle -f png

# IGV.js 页面(decoil_igv_config.py 生成配置 → track 渲染)
python modules/decoil/bin/decoil_igv_config.py \
    --sample s1:out/s1/decoil/reconstruct.ecDNA.filtered.bed \
    --bw s1:out/s1/s1.bigwig --igv-config igv_base.json \
    -o out/s1/visualization/igv_ecdna_config.json
python modules/track/bin/track.py --mode igv \
    --config out/s1/visualization/igv_ecdna_config.json \
    -i out/s1/s1.bigwig -o out/s1/visualization/igv_ecdna.html
```

## 2. 算法原理与生物学意义

ecDNA 是肿瘤中脱离染色体、自主复制的环状 DNA,携带高水平癌基因
扩增且结构高度异质。DeCoil(Giurgiu et al. 2024)利用长读长跨越
结构变异断点的证据,把 SV 组装成候选环:

- **reconstruct**:以 SV 断点为边、基因组片段为节点重建环/链,`summary.txt`
  每环一行给出拓扑(circular/linear 等)与占比(proportion);
  `--filter-score` 以上为高置信 ecDNA(`filtered.bed`);
- **pbsv2decoil.py**:pbsv VCF → sniffles1 字段(STRANDS/CHR2/END、
  FORMAT DR/DV),兼容上游 pbsv 调用;
- **decoil_plot**:环大小分布、拓扑构成、区段表(单样本 vs compare
  横向面板,含共有/独有区段);
- **decoil_structure**:逐环画环状结构(基因与覆盖叠放),无环时输出
  "No ecDNA detected" 占位图——阴性结果可交付;
- **igv/ucsc**:交互/公共浏览器轨道(bigwig + ecDNA BED + GTF),
  IGV 默认定位到最大 ecDNA 区段。

结果怎么读:`summary.txt` 的 proportion 是该环在样本 ecDNA 群体中的
占比;`reconstruct.ecDNA.filtered.bed` 的片段顺序即环的构成;结构图
上看断点两侧基因与扩增覆盖(高覆盖 = 高拷贝)。生物学意义:ecDNA
驱动癌基因超扩增与肿瘤内异质性,环构成与占比变化与耐药/复发演化
相关;跨样本共有/独有区段表回答「哪些环是共享的、哪些是私有的」。

## 3. 参考文献

1. Giurgiu I, et al. *Reconstructing extrachromosomal DNA structural
   heterogeneity from long-read sequencing data using Decoil*.
   Genome Research, 2024. DOI: 10.1101/gr.279123.124(方法出处;
   引文摘自模块 docstring)。
2. DeCoil —— 模块内未附仓库链接;可视化轨道脚本复用
   `modules/track/bin/track.py`(IGV.js / UCSC)。
