# exomePeak

MeRIP-seq/m6A-seq 甲基化峰调用与差异甲基化分析模块,封装 exomePeak R/Bioconductor 包
(2.17.0)。引擎脚本 `bin/exomePeak.r` 一条命令完成:峰调用/差异分析 → gene_id 翻译为
gene_name → 清理中间物,只交付注释后的 `*_gene_names.{bed,xls}`。

## 1. 模块输入输出与参数(怎么跑)

### 1.1 输入

| 输入 | 形态 | 说明 |
|---|---|---|
| IP/Input BAM | `{indir}/{sample_id}/{sample_id}{bam_suffix}.bam` | 形态 A(indir + 布局参数);默认 `bam_suffix=".dedup"`,即上游 igv 步 markdup 去重产物;call 与 diff 模式同一套输入构造(去重依据见 §3.5) |
| 样本分组 | config 列表 | ip/input 两组必填;treated_ip/treated_input 两组都非空才走差异分析,否则走峰调用 |
| 基因注释 | GTF 文件 | config `gtf`;取 `gene` 行的 gene_id/gene_name 建映射 |

### 1.2 参数(config 键 ↔ 引擎 CLI ↔ 默认值)

| config 键 | 引擎 CLI | 默认值 | 说明 |
|---|---|---|---|
| `indir` | — | `"input"` | 上游 BAM 目录(上游 config 的 `outdir`) |
| `outdir` | `--outdir` | `"output"` | 输出目录 |
| `logdir` | — | `"log"` | 规则日志目录,日志写 `{logdir}/exomePeak/<规则名>.log` |
| `ROOT_DIR` | — | `"."` | 项目根目录,用于定位 `bin/` 脚本 |
| `gtf` | `--gtf` | `""` | GTF 注释路径(必填) |
| `bam_suffix` | — | `".dedup"` | BAM 文件名布局后缀(`{sample_id}` 与 `.bam` 之间),call/diff 两规则共用同一构造 |
| `ip_samples` | `--ip_bams` | `[]` | 未处理组 IP 样本 ID 列表 |
| `input_samples` | `--input_bams` | `[]` | 未处理组 Input 样本 ID 列表 |
| `treated_ip_samples` | `--treated_ip_bams` | `[]` | 处理组 IP 样本 ID 列表 |
| `treated_input_samples` | `--treated_input_bams` | `[]` | 处理组 Input 样本 ID 列表 |
| `Procedure.Rscript` | — | `"Rscript"` | Rscript 可执行路径(可选覆盖) |

直接调用引擎(绕开 Snakemake)的示例:

```bash
Rscript bin/exomePeak.r \
    --gtf annotation.gtf \
    --ip_bams ip1.dedup.bam ip2.dedup.bam \
    --input_bams input1.dedup.bam input2.dedup.bam \
    --treated_ip_bams tip1.dedup.bam \
    --treated_input_bams tinput1.dedup.bam \
    --outdir results/
```

### 1.3 输出

`diff_exomePeak`(差异模式,处理组+未处理组):

| 文件 | 内容 |
|---|---|
| `diff_peaks_gene_names.{bed,xls}` | 全部检出峰 + 差异信息 |
| `sig_diff_peak_gene_names.{bed,xls}` | 显著差异甲基化峰 |
| `con_sig_diff_peak_gene_names.{bed,xls}` | 各重复一致的显著差异峰(**推荐交付集**) |

`call_exomePeak`(峰调用模式,单一条件):

| 文件 | 内容 |
|---|---|
| `all_peaks_gene_names.{bed,xls}` | 全部检出峰 |
| `con_peaks_gene_names.{bed,xls}` | 各重复一致峰(**推荐交付集**) |

中间物处理:exomePeak 原始表(`peak/diff_peak/sig_diff_peak/con_sig_diff_peak.*`,
映射见 §3.1)与 `exomePeak.Rdata` 落在 `{outdir}/exomePeak_output/`,全部转换成功后
该目录整体删除;`{outdir}/gene_map.rds` 是 GTF 解析缓存(签名 = 路径+大小+mtime),
保留以复用,GTF 变化自动失效。

命令回执:`{outdir}/diff_exomePeak.sh`(无时间戳,当前命令)+ `{OUTDIR_ROOT}/log/cmd/
diff_exomePeak_<时间戳>.sh`(历史备份);call 规则同理。

### 1.4 模块内脚本 CLI

`bin/geneId2name.r` 除被引擎引用外可独立调用:

| 旗标 | 必填 | 说明 |
|---|---|---|
| `-i/--infile` | 是 | 输入表(bed/xls,需含 `name` 列) |
| `--gtf` | 是 | GTF 注释 |
| `-o/--outfile` | 是 | 输出表(追加 `gene_name` 列) |
| `--gene_id_col` | 否 | gene_id 所在列名,默认 `name` |

## 2. 结果解读(怎么读)

BED 与 XLS 前 12 列为标准 BED12(https://genome.ucsc.edu/FAQ/FAQformat.html#format1):
chrom、chromStart、chromEnd、name(gene_id)、score(峰 p 值)、strand、thickStart、
thickEnd、itemRgb、blockCount、blockSizes、blockStarts;峰可跨多个外显子。

XLS 在 BED12 之后追加:

| 列 | 含义 |
|---|---|
| `lg.p` | log10(p 值),峰作为甲基化位点的显著性 |
| `lg.fdr` | log10(FDR),多重检验校正后显著性 |
| `fold_enrchment` | 峰区内 IP/Input 富集倍数 |
| `diff.lg.fdr` | log10(FDR),差异甲基化显著性(仅 diff 输出) |
| `diff.lg.p` | log10(p 值),差异甲基化显著性(仅 diff 输出) |
| `diff.log2.fc` | log2(odds ratio);>0 处理组高甲基化,<0 低甲基化(仅 diff 输出) |
| `gene_name` | 由 gene_id 翻译的基因名,逗号分隔多个;翻译不到的保留原 gene_id |

选集建议:重复间一致集(`con_*`)可重复性最高,diff 模式交付 `con_sig_diff_peak_gene_names.*`。

## 3. 算法原理(怎么算)

### 3.1 峰调用

按 GTF 注释把转录组切分为 batch,在窗口(WINDOW_WIDTH=200 nt,SLIDING_STEP=30 nt)
内统计 IP/Input reads 计数,对每个窗口做单侧二项检验(IP 计数在 IP+Input 总计数下
偏离 Input 背景比例的概率),合并相邻显著窗口成峰(最短 50 nt),要求
FDR < 0.05 且 fold enrichment > 1;reads 比对质量过滤 MAPQ ≥ 30,并剔除局部
标签异常(REMOVE_LOCAL_TAG_ANOMALITIES)。原始表写入 `{outdir}/exomePeak_output/`:
call 模式为 `peak.{bed,xls}` 与 `con_peak.{bed,xls}`,diff 模式为 `diff_peak.*`、
`sig_diff_peak.*`、`con_sig_diff_peak.*`;引擎把它们逐一翻译为 `*_gene_names.*`
交付表后连同 `exomePeak.Rdata` 删除。注意:双端数据按单端处理,配对信息被忽略 [1]。

### 3.2 一致峰(consistent peaks)

对每个重复单独评价富集,各重复均满足 p < 0.05 且 fold enrichment > 1 的峰进入一致集
(`con_*`),表示在全部 IP 重复中稳定富集。

### 3.3 差异甲基化分析

四组文库(IP/Input × 处理/未处理)计数先按文库大小相互归一,再对 2×2 计数表做
超几何检验(双侧取较小尾概率),`diff.log2.fc` = log2((处理 IP/处理 Input)/(未处理
IP/未处理 Input)) 即 log2 odds ratio;BH 法 FDR 校正(总 reads ≤ 10 的位点不参与
FDR 估计)。显著差异:diff FDR < 0.05 且 |fold change| > 1;一致差异再叠加各重复
单独方向一致(p < 0.05、|FC| > 1、符号一致)。备选检验:卡方检验(`DIFF_PEAK_METHOD`)。

### 3.4 gene_id → gene_name

`bin/geneId2name.r` 流式解析 GTF 的 `gene` 行建 gene_id→gene_name 映射,读入峰表将
`name` 列逐个(逗号分隔多个)翻译,追加 `gene_name` 列;映射缓存为 `gene_map.rds`。

### 3.5 为什么用去重 BAM

MeRIP-seq 峰调用的标准流程为 QC → 比对 → 过滤(去 rRNA 与 PCR 重复)→ 峰调用,
PCR 重复去除是过滤步的常规组成 [3];应用去重的流程(如 CVm6A 用 Picard 去重)峰集
更严格、假阳性更少 [3]。exomePeak 官方协议 [1] 未要求保留重复 reads。本模块 call 与
diff 模式统一消费上游 igv 步 markdup 产物 `{sample_id}.dedup.bam`,保证两种分析
串度一致,且该文件是上游规则声明的 output(DAG 依赖可追踪)。

## 4. 参考文献

1. Meng J, Lu Z, Liu H, Zhang L, Zhang S, Chen Y, Rao MK, Huang Y. A protocol for RNA methylation differential analysis with MeRIP-Seq data and exomePeak R/Bioconductor package. Methods. 2014;69(3):274-281. doi:10.1016/j.ymeth.2014.06.008. PMID: 24979058.
2. Zhou J, Wei Z, Zhen D, Wang Y, Su J, Coenen F, Jia G, Rigden DJ, Meng J. Comprehensive Epitranscriptome Analysis from MeRIP-seq Data with exomePeak2. Genomics Proteomics Bioinformatics. 2026. doi:10.1093/gpbjnl/qzag019. PMID: 41776836.
3. Song R, Sutton GJ, Li F, Liu Q, Wong JJL. Variable calling of m6A and associated features in databases: a guide for end-users. Brief Bioinform. 2024;25(5):bbae434. doi:10.1093/bib/bbae434. PMID: 39258883.
4. exomePeak R/Bioconductor package (Jia Meng, 2019). https://github.com/ZW-xjtlu/exomePeak
