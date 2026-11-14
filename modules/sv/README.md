# sv — 结构变异下游分析模块(比较特异、MEI、circos、OncoPrint 等)

对 pbsv 产出的 SV VCF 做队列/分组下游分析:每对比提取比较特异 SV、
VEP + RepeatMasker 注释、统计与 TE 富集,MEI(移动元件插入)指纹分析,
circos 图、基因模型图、跨组差异比较与 OncoPrint。方法文档(输入契约、
判定规则、参数依据、方法学教训)随模块维护在 `docs/`(MEI_METHODS.md、
README_polyA_TSD.md、TRA_TE_ANALYSIS_METHODS.md),本 README 只做
接线说明;polygenomes 变体见 `polygenomes/`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `comparisons` (config,顶层) | DESeq2 风格对比定义字典:`{"{control}_vs_{experimental}": {control_group_name, experimental_group_name, control_samples: [...], experimental_samples: [...]}}`;wildcard `{comparison}` 取其键 |
| `vcf_dir` + `vcf_suffix` (config) | 每样本 SV VCF(pbsv 布局):`{vcf_dir}/{sample}/{sample}{vcf_suffix}`(默认 `.sv.vcf.gz`);路径不进 config |
| `genome.fasta` / `genome.gtf` (config) | 参考 FASTA 与基因注释(**必填**,`_require` 校验) |
| `genome.custom_gtf` / `custom_fasta` / `te_lib` / `repeat_bed` (config) | 可选:自定义注释、TE 共识库(legacy 分类器)、RepeatMasker BED(sv_mei) |
| `Params.sv.mei.dfam_h5` (config) | Dfam FamDB h5 分区(列表,MEI 分类器首选;与 `genome.te_lib` 二者须有其一) |
| `Params.sv.oncoprint.cosmic_file` (config) | COSMIC 癌基因列表(必需,`_require`) |
| `Params.sv.oncoprint.deseq2_files` (config) | `{condition: path}` DESeq2 结果字典(OncoPrint 叠加表达) |

### 1.2 输出

```
# 每 comparison(每对比)
{outdir}/{comp}/{comp}_merged.vcf            # SURVIVOR 合并 VCF(sv_exp_specific)
{outdir}/{comp}/{comp}_only.vcf              # 比较特异 SV
{outdir}/{comp}/{comp}_annotated.tab         # VEP 注释表
{outdir}/{comp}/repeatmasker/                # RepeatMasker 注释产物
{outdir}/{comp}/statistics/                  # 统计表
{outdir}/{comp}/enrichment/                  # TE 富集分析
{outdir}/{comp}/circos/                      # circos 图(sv_exp_circos)
{outdir}/{comp}/gene_model/                  # 基因模型 + SV 叠加图(sv_gene_model)
{outdir}/{comp}/mei/mei.done                 # MEI 分析完成标记(产物与缓存在 mei/ 内)

# 跨对比
{outdir}/sv_diff_analysis/                   # 组间 SV 差异分析(pbsv_sv_diff_analysis.py)
{outdir}/sv_oncoprint/oncoprint_matrix.csv   # OncoPrint 矩阵(+ sv_oncoprint/ 图)
```

回执:`{outdir}/{comp}/sv_*.sh`、`{outdir}/sv_{diff_analysis,oncoprint}.sh`;
日志 `{logdir}/sv/{comp}/sv_*.log` 等。`sv_result` 聚合全部产出
(enrichment 目录由 sv_exp_specific 一并生成)。

### 1.3 参数

config 位于 `Params.sv.<段>` 下(**注意**:`dist` 例外,读
`Params.exp_specific.dist`;示例 json 把各段放在 `Params.<段>` 而无 `sv`
层——以 .smk 读取路径为准):

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.image_formats` | `-f`(可重复) | `[png]` | 全部出图格式 |
| `Params.genes` | `-g`(gene_model,可重复) | [] | 基因模型图的目标基因 |
| `Params.sv.exp_specific.vep_cache` | `--vep_cache` | `~/.vep` | VEP 缓存目录 |
| `Params.sv.exp_specific.species` / `.assembly` | `--species` / `--assembly` | `mus_musculus` / `GRCm39` | VEP 物种/组装 |
| `Params.exp_specific.dist` | `-d` | 500 | 特异 SV 判定距离(bp) |
| `genome.custom_gtf`(+`.custom_fasta`) | `--custom_gtf` / `--custom_fasta` | "" | 非空才传 |
| `genome.te_lib` | `--te_lib` / `--te-lib` | "" | TE 共识库(exp_specific 富集 / mei legacy 分类器) |
| `Params.sv.circos.genome` / `.cytoband` / `.ins_bin_size` | `--genome` / `--cytoband` / `--ins_bin_size` | `mm39` / "" / 100000 | circos 组装名、染色体带、INS 分箱 |
| `Params.sv.gene_model.threads` | `-j` | 1 | 并行画图 |
| `Params.sv.diff_analysis.large_sv_threshold` | `-s` | 10000 | 大 SV 长度阈值(bp) |
| `Params.sv.oncoprint.cosmic_file` | `--cosmic_file` | — | 必需 |
| `Params.sv.oncoprint.deseq2_files` | `-d cond:path` | {} | 表达叠加 |
| `Params.sv.mei.min_len` | `--min-len` | 100 | MEI 最小 INS 长度 |
| `Params.sv.mei.full_length` / `.full_length_family` | `--full-length-l1` / `--full-length-family` | 6000 / `L1` | 全长 L1 判定阈值/家族 |
| `Params.sv.mei.dfam_h5` | `--dfam-h5`(可重复) | [] | Dfam FamDB h5(首选分类器) |
| `Params.sv.mei.rm_species` / `.dfam_mode` | `--rm-species` / `--dfam-mode` | "" / `export` | RepeatMasker 物种/FamDB 导出模式 |
| `Params.sv.mei.te_buckets` | `--te-buckets`(JSON) | 见示例 json | TE 分桶定义(Alu/L1/SVA/LTR/DNA/microsatellite/other/unclassified) |
| `Params.sv.mei.force_rerun` / `.enabled` | `--force-rerun` / — | false / true | 强制重算;`enabled` 供 node.py 注册门控 |
| (threads) | — | 4(sv_mei) 等 | `rule_res` 可覆盖 |
| (mem_mb) | — | 8000 / 16000(mei) | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致,以 `young_vs_aged` 为例):

```bash
# 比较特异 SV:合并 → 提取 → VEP/RepeatMasker 注释 → 统计/TE 富集
python modules/sv/bin/exp_specific.py \
    -c young1.sv.vcf.gz young2.sv.vcf.gz \
    -e aged1.sv.vcf.gz aged2.sv.vcf.gz \
    -o out/young_vs_aged/young_vs_aged \
    -d 500 --vep_cache ~/.vep --species mus_musculus --assembly GRCm39 \
    --annotate_format tab -f png

# MEI 分析(Dfam 分类器)
python modules/sv/bin/mei_analysis.py \
    -i out/young_vs_aged/young_vs_aged_merged.vcf \
    -o out/young_vs_aged/mei \
    --min-len 100 --full-length-l1 6000 --full-length-family L1 \
    --threads 4 --title young_vs_aged \
    --dfam-h5 /path/to/dfam_famdb.h5 --rm-species "Mus musculus" \
    --dfam-mode export --te-buckets '[...]' \
    --ctrl-name young --exp-name aged

# 跨组差异 + OncoPrint
python modules/sv/bin/pbsv_sv_diff_analysis.py -o out/sv_diff_analysis \
    -s 10000 -g young_vs_aged:out/young_vs_aged/young_vs_aged_only.vcf -f png
python modules/sv/bin/run_OncoPrint.py --cosmic_file cosmic_cancer_genes.txt \
    -o out/sv_oncoprint -g young_vs_aged:out/young_vs_aged/young_vs_aged_annotated.tab \
    -d ctrl:deseq2_ctrl.tsv -f png
```

## 2. 算法原理与生物学意义

各分析线(`bin/` 脚本,方法细节见 `docs/`):

- **sv_exp_specific**:SURVIVOR 合并组内 VCF → 提取仅存于实验组(或
  对照组)的比较特异 SV(`-d` 距离内判同)→ VEP 功能注释 +
  RepeatMasker 重复注释 → 统计(SV 类型/大小分布)→ TE 富集(断点
  是否富集于 TE);
- **sv_mei**:对合并 VCF 的 INS 等位基因按 RepeatMasker/Dfam 分类
  (te_buckets 分桶:Alu/L1/SVA/LTR/DNA/微卫星/其他),检测 poly(A)
  与 TSD(target site duplication)逆转录转座指纹(判定逻辑与 TPRT
  生物学原理见 `docs/README_polyA_TSD.md`),标注基因组上下文并出
  6 图套件;TE 插入年代测定(te_age_analysis.py)是独立分析,不在此
  运行;
- **sv_exp_circos / sv_gene_model**:特异 SV 的全基因组 circos 视图
  (INS 分箱)与指定基因(`Params.genes`)的基因结构 + SV 断点叠加图;
- **sv_diff_analysis**:`-s` 阈值划分大/小 SV,跨对比比较 SV 负荷与
  类型构成;
- **sv_oncoprint**:按 COSMIC 癌基因列表汇总各组 SV 命中,叠 DESeq2
  表达,输出 oncoprint 矩阵与图——驱动事件的队列级总览。

生物学意义:长读长 SV 揭示重复区结构变异;MEI/TRA 等分析把断点与
TE 联系起来(TE 介导的重排、逆转录转座插入),OncoPrint 汇总驱动
基因命中。证据分级(T0-T3)与局限见 docs 约定:不确定结论只标
「候选/假设」。

## 3. 参考文献

1. 方法文档(与代码同版本):`docs/MEI_METHODS.md`、
   `docs/README_polyA_TSD.md`、`docs/TRA_TE_ANALYSIS_METHODS.md`
   (含参数依据与方法学教训)。
2. Ensembl VEP —— https://www.ensembl.org/info/docs/tools/vep/index.html
3. RepeatMasker —— https://www.repeatmasker.org;Dfam —— https://www.dfam.org
4. COSMIC —— https://cancer.sanger.ac.uk/cosmic(oncoprint 基因列表来源)
5. 运行环境:`sv.yaml` / `sv.def`(sv.sif:RepeatMasker/blastn/pysam/
   pandas/matplotlib;samtools 用独立 SIF,见 `docs/README.md`)。
