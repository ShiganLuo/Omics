# sv — 结构变异下游分析模块(比较特异、MEI、circos、OncoPrint 等)

对 pbsv 产出的 SV VCF 做队列/分组下游分析:每对比提取比较特异 SV、
VEP + RepeatMasker 注释、统计与 TE 富集,MEI(移动元件插入)指纹分析,
circos 图、基因模型图、跨组差异比较与 OncoPrint。方法文档(输入契约、
判定规则、参数依据、方法学教训)收录于本 README §3(MEI 管线 §3.1、
poly(A)/TSD 指纹 §3.2、TRA/TE-TE 重组 §3.3、年轻插入验证 §3.4),§1-§2 为接线说明;
polygenomes 变体见 `polygenomes/`。

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
{outdir}/{comp}/tra_te/tra_te.done          # TE 介导易位分析完成标记(5 TSV + figs/ 4 图)

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
| `Params.sv.formats` | `-f`(可重复) | `[png]` | 全部出图格式(matplotlib 全集;旧键 `Params.sv.image_formats` 兼容读取,以 `formats` 为准) |
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
| `Params.sv.tra_te.n_perm` | `--n-perm` | 1000 | TRA 断点 TE 富集的置换次数(单断点采样;双断点对数固定 2×n_perm) |
| `Params.sv.tra_te.seed` | `--seed` | 42 | 置换随机种子(可复现) |
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
    --ctrl-name young --exp-name aged -f png

# TE 介导易位(TRA)分析(T0-T3 漏斗, 出 5 TSV + 4 图)
python modules/sv/bin/utils/tra_te_analysis.py \
    --vcf out/young_vs_aged/young_vs_aged_merged.vcf \
    --te-bed ref/rmsk_TE.lifted.bed \
    --fasta ref/genome.fa \
    --fai ref/genome.fa.fai \
    -o out/young_vs_aged/tra_te \
    --threads 4 --n-perm 1000 --seed 42

# 跨组差异 + OncoPrint
python modules/sv/bin/pbsv_sv_diff_analysis.py -o out/sv_diff_analysis \
    -s 10000 -g young_vs_aged:out/young_vs_aged/young_vs_aged_only.vcf -f png
python modules/sv/bin/run_OncoPrint.py --cosmic_file cosmic_cancer_genes.txt \
    -o out/sv_oncoprint -g young_vs_aged:out/young_vs_aged/young_vs_aged_annotated.tab \
    -d ctrl:deseq2_ctrl.tsv -f png
```

## 2. 算法原理与生物学意义

各分析线(`bin/` 脚本,方法细节见 §3):

- **sv_exp_specific**:SURVIVOR 合并组内 VCF → 提取仅存于实验组(或
  对照组)的比较特异 SV(`-d` 距离内判同)→ VEP 功能注释 +
  RepeatMasker 重复注释 → 统计(SV 类型/大小分布)→ TE 富集(断点
  是否富集于 TE);
- **sv_mei**:对合并 VCF 的 INS 等位基因按 RepeatMasker/Dfam 分类
  (te_buckets 分桶:Alu/L1/SVA/LTR/DNA/微卫星/其他),检测 poly(A)
  与 TSD(target site duplication)逆转录转座指纹(判定逻辑与 TPRT
  生物学原理见 §3.2),标注基因组上下文并出
  6 图套件;TE 插入年代测定(te_age_analysis.py)是独立分析,不在此
  运行;
- **sv_tra_te**:TE 介导易位(TRA)分析,T0-T3 证据漏斗(方法见
  §3.3)——合并 VCF 提取跨染色体 TRA 并镜像去重为独立事件(T0),
  断点落 RepeatMasker TE 区间判定 TE 相关(T1),置换零分布 + 二项
  检验验证双断点 TE 富集、官方 family_id 配对 vs shuffle 期望(T2),
  两断点 TE 拷贝序列 blastn 同源验证、NAHR 级判定(ident≥90% 且
  HSP(High-scoring Segment Pair,高分配对片段,blastn 局部比对块,
  §7.1)≥100bp,T3),附 4 图套件(enrichment/family pairs/homology/
  panorama);
- **sv_exp_circos / sv_gene_model**:特异 SV 的全基因组 circos 视图
  (INS 分箱)与指定基因(`Params.genes`)的基因结构 + SV 断点叠加图;
- **sv_diff_analysis**:`-s` 阈值划分大/小 SV,跨对比比较 SV 负荷与
  类型构成;
- **sv_oncoprint**:按 COSMIC 癌基因列表汇总各组 SV 命中,叠 DESeq2
  表达,输出 oncoprint 矩阵与图——驱动事件的队列级总览。

生物学意义:长读长 SV 揭示重复区结构变异;MEI/TRA 等分析把断点与
TE 联系起来(TE 介导的重排、逆转录转座插入),OncoPrint 汇总驱动
基因命中。证据分级(T0-T3)与局限见 §3.0 约定:不确定结论只标「候选/假设」。

## 3. 方法文档(与代码同版本, 供审核/复现)

只记录方法: 输入契约、算法、判定规则、统计设计、参数依据、质控点、
方法学教训、局限。**不含具体数据结果**(结果在数据输出目录)与
**运行用法**(用法见代码 argparse 帮助与 sv.smk 规则调用)。

| 节 | 覆盖分析 | 对应代码入口 |
|---|---|---|
| §3.1 MEI 管线方法 | 输入契约/五步流程/判定逻辑三层/差异分析/局限 | `bin/mei_analysis.py`(**唯一入口**, 分析→图一条龙; 绘图经 `bin/utils/mei_plot.py` 自动串联) |
| §3.2 poly(A)/TSD 指纹 | 指纹检测算法 + TPRT 生物学原理 | `bin/mei_analysis.py` 的 fingerprint 步骤 |
| §3.3 TRA/TE-TE 重组 | T0-T3 证据漏斗、统计设计、SURVIVOR 行为定性、序列级同源验证(含序列背景规格) | `bin/utils/tra_te_analysis.py` |
| §3.4 年轻插入验证 | INS 转座判定、定年筛选、VAF 验证 | `bin/utils/te_age_analysis.py` + `bin/mei_analysis.py` 事件表 |

### 3.0 约定

- 结论标注证据级(T0-T3 等), 不确定的写"候选/假设", 不下硬结论;
- 参数给出默认值与依据;
- 历次方法学错误的修正记录保留在文档内("方法学教训"节), 供审核追溯。

### 3.1 MEI (Mobile Element Insertion) 分析方法文档

针对 SURVIVOR 合并 SV VCF(pbsv 来源)的转座事件分析管线的方法说明。
只记录方法(输入契约、算法、判定规则、参数依据、局限), 不含具体数据结果
与运行用法; 运行方式见代码入口的 argparse 帮助与工作流调用
(`modules/sv/sv.smk` sv_mei 规则)。
代码入口: `modules/sv/bin/mei_analysis.py`(唯一入口; 绘图经
`bin/utils/mei_plot.py` 串联)。TE 插入定年(divergence-from-consensus
分子钟)为独立分析(`bin/utils/te_age_analysis.py`), 不在本管线内。
poly(A)/TSD 指纹的检测算法与 TPRT 生物学原理详见 §3.2。

#### 一、输入契约

| 输入 | 要求 | 方法学作用 |
| --- | --- | --- |
| merged VCF | SURVIVOR merge(pbsv 单样本 calls) | INS 的 ALT 携带完整插入序列(pbsv AAL 字段被 SURVIVOR 保留), 是 TE 分类的原始证据 |
| TE 一致序列 FASTA | 物种 TE 库 | RepeatMasker 分类的比对库 |
| 基因注释 GTF(.gz) | gene/exon 层 | exon/intron/intergenic 语境注释 |
| TE 区域 BED | rmsk lift-over | TE-into-TE 判定 |
| 参考 FASTA(已索引) | 需 .fai | TSD 检测取断点侧翼 |

SUPP_VEC 语义: 对照组在前、实验组在后的 merge 顺序下, `10`=对照特异、
`01`=实验特异、`11`=共享; 多样本时逐位对应。读取时 `01` 会丢前导零变
`1`, 必须补零对齐(此为已踩过的解析坑)。

#### 二、分析流程(五步)

1. **INS 提取**: 解析 VCF, 取 SVTYPE=INS 且 ALT 带序列的记录(REF/ALT
   共享前缀裁剪后得纯插入序列); symbolic `<INS>` 无序列者丢弃并计数。
2. **TE 分类**: 插入序列 ≥阈值(默认 100bp)写 FASTA → RepeatMasker
   `-s -no_is -lib <te_lib>`; 每条取**覆盖跨度最大**的命中(注意不是 SW
   最高), 记录亚家族/覆盖度/链方向。家族归桶: Alu(SINE)、L1(LINE)、
   SVA、LTR、DNA、simple_repeat、other、unclassified。低于阈值者标
   short_unclassified(多为测序噪声, 保留不删)。.out 解析用"右锚定
   回退行走"(从行尾数字块回溯定位 class/repeat/strand 字段)防止列错位。
3. **MEI 指纹**: poly(A) 尾与 TSD 检测, 算法与原理见 `README_polyA_TSD.md`。
   判读不对称性: 阳性可确证 TPRT 插入, 阴性不能证伪。
4. **语境注释**: GTF gene/exon → exon/intron/intergenic; rmsk BED →
   是否落入已有 TE 区域(TE-into-TE vs 独特位点)。
5. **全长 L1 候选**: SVLEN≥阈值(默认 6000)且家族=L1(潜在有活性元件
   清单; 未做 ORF 完整性检查, 需要时可 ORFfinder 交叉验证)。

#### 三、判定逻辑(三层区分, 勿混用)

1. **序列来源层**(直接观测): 插入序列与哪个 TE 家族一致序列同源
   (RepeatMasker best-hit)。只能得出"序列是 TE 来源/重复来源/未知",
   **不能直接等同于转座事件**。
2. **移动机制层**(文献依据, 按家族不同):
   - Alu / L1 / SVA: L1 ORF2p 介导的 TPRT, 复制型(RNA 中间体)——
     polyA/TSD 指纹只对此机制有效;
   - LTR / ERV: 整合酶介导的类逆转录病毒机制, 复制型但不经 TPRT,
     无 polyA/TSD 指纹;
   - DNA 转座子: 转座酶介导的剪切-粘贴(cut-and-paste), 非复制型;
     哺乳动物基因组中多为化石拷贝(无活性)。
   - simple_repeat / microsatellite: 复制滑动/不等交换的串连重复扩张,
     不是转座。
3. **事件层**(推断, 需证据分级):
   - 强: TE 来源 + polyA/TSD 指纹 + 低分歧度 → TPRT 新生插入(转座事件);
   - 中: TE 来源 + 群体多态 → 历史转座事件的遗传残留(产物在, 事件古老);
   - 弱: 仅 TE 同源性, 不排除重复区比对伪影。

**未分类的处理(声明为假设, 非结论)**: other/unclassified 性质未知。
若假设均非 TE, 得转座占比下界; 若全为 TE, 得上界。一切"转座占比"
结论按区间理解, 不给单点值。

#### 四、输出 schema(接口定义, 非用法)

per-event 表(坐标/长度/组别/家族/亚家族/覆盖度/polyA/TSD/语境/基因/
TE 区域/全长 L1 标记)、家族×组别计数矩阵、语境计数、全长 L1 清单、
汇总 json、图组(family_by_group / insert_length_by_family /
genomic_context / te_region_ratio / full_length_L1 /
polya_tsd_fingerprint)、自包含 HTML 报告、RepeatMasker 原始输出留档
(可复查亚家族)。具体文件名见代码 docstring。

#### 五、方法学局限

1. **1v1 设计**: 每组单样本时"组特异"是个体多态而非组效应; 组间定量
   比较需每组多重复(SUPP_VEC 逐位自适应, 方法不变)。
2. **polyA 检出率是下限**: caller 对 ALT 序列有截断, 3' 端常缺失,
   人为压低检出率; 长元件(如 L1, 多为 5' 截短拷贝)受影响更大。
3. **TSD 仅阳性可确证**: 参考位点已有祖先插入时断点两侧无重复序列;
   且断点 1bp 偏移即配对失败, 真实 SV call 上检出率天然极低——
   TSD 面板是"能确证的下限", polyA 才是主判据。
4. **RepeatMasker 分类局限**: 比对的是共识序列而非基因组拷贝; 嵌合/
   重排插入只报最长命中; 低覆盖度命中分类可信度低。

#### 六、差异分析(mei_diff, 探索性)

在 ctrl/exp 两分支间给 TE 家族/亚家族排序, 回答"哪个最不一样"。无
生物学重复, 结论仅限探索性(见五.1)。

**分组单位**: `mei_diff.tsv`(family 级 = TE 桶)与
`mei_diff_by_repeat.tsv`(te_repeat 级 = RepeatMasker 亚家族名, 无命中
记录不参与)。summary.json 另存 `diff_top` / `diff_top_repeat`(前 10)。

**表列语义**:

- `n_ctrl_specific` / `n_exp_specific` / `n_shared`: 两分支特异调用数与
  共享数(shared 计入两分支的指标计算);
- `log2fc_exp_vs_ctrl`: 火山图横轴,
  `log2(((b+0.5)/(N_exp+1)) / ((a+0.5)/(N_ctrl+1)))`——该家族占 exp 特异
  调用的**份额**与占 ctrl 特异调用份额之比的对数(a/b = 该家族两分支特异
  数, N = 各分支特异总数)。用份额而非原始计数, 以约掉两分支总调用数的
  不对称; 伪计数防小计数爆炸; x = ±1 即份额 2 倍差;
- `p_value`: 双侧 Fisher 精确检验, 2×2 表 `[[a, N_ctrl-a], [b, N_exp-b]]`。
  零假设 = 该家族特异调用按两分支总量比例随机分配; **边缘固定**使全局
  调用率差异被条件化进零假设, 且稀疏格(4 vs 0)下依然精确(卡方不适用);
- `fdr`: BH 校正(每张表内多重检验);
- `ctrl_/exp_polyA_pct`、`in_te_pct`、`median_ins_len`: 两分支的 poly(A)
  率、落入参考 TE 区率、中位插入长度。

**重点关注规则(`focus` 列, 候选筛选非显著性声明)**: 同时满足

1. `|log2FC| >= 1`——占比 2 倍差(效应量门槛);
2. `p < 0.05`——失衡超出随机分配(水平线, 图中已标);
3. `n_ctrl_specific + n_exp_specific >= 10`——排除小计数噪声比值。

**火山图(mei_diff_volcano.png)**: x = log2FC, y = -log10(p), 不画参考
线; 显著点(p<0.05)黑边、标签加 `*`; 标注 = 负数区间(ctrl 富集)与正数
区间(exp 富集)各前 5 个, 按联合分数 `|log2FC| * -log10(p)` 排序("双高"
规则——单维极端(大计数小效应、或大效应弱证据)不入选), 标签自动避障。
图只画 te_repeat 级; 家族级查 mei_diff.tsv。

**判读要点**: p 随计数放大——大计数家族的微小失衡也会"显著"(如
ALR/Alpha 型技术失衡), 先看 |log2FC| 与计数规模, 再看 p; focus 候选需
用 rm_id 回查多命中结构(完整前病毒 vs 片段)、polyA/in_te 列, 并以独立
证据(RNA-seq、重复样本)复核后才能作为生物学结论。

### 3.2 MEI 指纹检测原理:poly(A) 尾与 TSD —— 为什么它们是 L1 介导插入的标记

只记录原理、算法与判定依据; 不含数据结果与运行用法。
检测实现位于 `modules/sv/bin/mei_analysis.py`, 管线总览见 §3.1。

---

#### 一、背景:MEI 与 L1 介导的逆转座

MEI(Mobile Element Insertion,移动元件插入)指转座元件通过"复制-粘贴"
机制在基因组中产生的新插入。人类/灵长类基因组中活跃的 MEI 家族有三个:

| 家族 | 类型 | 大小 | 自身是否编码逆转录酶 |
| --- | --- | --- | --- |
| L1 | LINE | ~6 kb,全长含 ORF1/ORF2 | 是(ORF2p = 内切酶 + 逆转录酶) |
| Alu | SINE | ~300 bp | 否,完全依赖 L1 机器 |
| SVA | SVE | ~0.7–4 kb | 否,完全依赖 L1 机器 |

关键点:**Alu 和 SVA 是非自主元件**。它们自身只有 RNA(借宿主 Pol III /
Pol II 转录),没有编码任何移动所需的蛋白,必须"劫持"细胞中 L1 表达的
ORF2p 才能完成插入(in trans mobilization)。因此:

- L1 的插入由自己的 ORF1p/ORF2p 介导(顺式偏好,但也能反式动员别人);
- Alu/SVA 的插入由 **L1 的 ORF2p** 介导。

三者用的是同一部机器、同一条插入通路,所以在插入位点留下同一套结构指纹。
这就是"poly(A) 尾 + TSD = L1 介导插入的标记"这句话的完整含义——它不仅
标记 L1 自身的插入,也标记被 L1 机器反式动员的 Alu/SVA 插入。

对比:LTR 逆转座子(类逆转录病毒,经 RNA 中间体但走 cDNA 剪接/整合酶
途径)和 DNA 转座子(剪切-粘贴)不走这条通路,不产生这套指纹。这也是
检测只对 family ∈ {L1, Alu, SVA} 做 TSD 判定、指纹图只画这三个家族的原因。

#### 二、TPRT:指纹的产生机制

L1 介导的插入遵循 **TPRT(Target-Primed Reverse Transcription,靶标引物
逆转录)** 机制(Feng et al. 1996, Cell;Luan et al. 1993, Cell),分步如下:

```
1. 转录        基因组中源拷贝 → L1/Alu/SVA RNA,3' 端带 poly(A) 尾
                  (来自转录终止时 Pol 的滑动合成,或 A-rich 3' 尾)
2. 切口        ORF2p 的内切酶结构域在基因组 T 富集位点切割:
                  5'-...TTTT/AAAA...-3'  偏好切在 T(A) 丰富区,
                  切口产生的 3'-OH 暴露
3. 引物退火    RNA 的 poly(A) 尾与切口处 T 富集序列互补配对
4. 逆转录      ORF2p 以该 3'-OH 为引物、RNA 为模板,在基因组上
                  直接合成 cDNA(第一链);"copy-and-paste"的核心步骤
5. 第二链      互补链合成、连接、修复
6. 结果        插入序列 = 元件 cDNA + 3' 端 poly(A) 尾
                  断点两侧 = 同一段基因组序列各留一份(TSD)
```

两个指纹的来源:

- **poly(A) 尾**:第 1/3/4 步——RNA 上的 poly(A) 尾既是退火"把手",
  也被原样逆转录成 cDNA,固化在插入序列的 3' 端。
- **TSD(Target Site Duplication,靶位点重复)**:第 2/6 步——ORF2p
  切口时,断点上游那一小段基因组序列(典型 7–20 bp,长度与元件/物种
  相关,L1 约 15–20 bp,Alu 约 7–15 bp)在修复完成后被复制到插入位点
  的两侧,形成首尾相同的直接重复:

```
插入前:  5' —————— ABCDEFGH —————— 3'          (ABCDEFGH 为基因组原有序列)
插入后:  5' ——— ABCDEFGH [元件 + polyA] ABCDEFGH ——— 3'
                 ↑                                ↑
              5' TSD                          3' TSD
```

另外 TSD 两侧常伴随靶位点**微缺失/微同源**,插入常带 5' 截短(L1 尤甚,
因为逆转录经常未完成),这些共同构成 L1 插入位点的"化石记录"。

判读时的不对称性(重要):

- **阳性即可确证**:插入序列末端检出 A-rich 尾、断点两侧检出相同 k-mer
  重复,几乎只能用 TPRT 解释,是 L1 介导新生插入的直接证据。
- **阴性不能证伪**:polyA 可能因 caller 对 ALT 序列截断而丢失;参考基因
  组该位点若已带祖先插入(群体/物种固定拷贝),则参考中根本不存在"插入
  前"的侧翼序列,TSD 必然检不出。所以 TSD/polyA 缺失 ≠ 非转座事件。
  指纹检出率应作为"携带指纹的插入占比",即 L1 介导插入比例的**下限**。

#### 三、检测算法与判定依据

##### 3.1 poly(A) 尾检测

- 窗口:插入序列**末 60 bp**(不足 60 取全长)。poly(A) 尾长度典型
  20–150 bp,60 bp 窗口兼顾灵敏度与特异性(全序列统计会被长元件稀释)。
- 判阳条件(满足其一):
  1. 窗口内**连续 A ≥12 bp**——真实测序/比对错误率下,随机序列出现
     12 连 A 的概率 ≈ (1/4)^12 ≈ 6e-8,基本排除背景;
  2. 窗口内 **A 占比 ≥90%**——容忍测序错误造成的零星打断,覆盖
     "A-rich 但不纯"的真实尾巴(逆转录/复制过程中的滑移、替换)。
- 前置条件:插入长度 ≥200bp 才检测(太短的插入多为噪声,记缺失)。
- 只查 3' 端、不区分链方向:MEI 的 ALT 序列方向由 caller 决定,polyA
  若在序列 5' 端说明该序列是反向的,靠 RepeatMasker 的链方向列交叉判断。

##### 3.2 TSD 检测

- 坐标约定:VCF 中 INS 的 POS 是插入前 1-based 锚定碱基,插入发生在
  POS 与 POS+1 之间;pysam fetch 是 0-based 半开区间,故上游取
  `[pos-1-flank, pos)`,下游取 `[pos, pos+flank+1)`(flank=30)。
- 判定标准:紧邻断点上游末 k bp 与下游首 k bp**完全相同**(k 从 20
  递减到 7)。要求"紧邻 + 两侧同时存在 + 精确匹配",对应 TPRT 产生的
  典型 7–20 bp 直接重复;全 N(参考缺口)排除。
- 从长到短扫描:一个真实 TSD 若两侧完好,最长匹配就是它本身;先扫短
  k 会把 TSD 的子串误当 TSD,长度信息也丢失。
- 只对 family ∈ {L1, Alu, SVA} 检测(LTR/DNA 转座子无此机制,检测无
  意义,且会引入假阳性)。
- 输出语义:重复序列本身与长度分列记录,未检出为空。

##### 3.3 指纹率统计

polyA 率 = family 内 polyA 阳性比例;TSD 率 = family 内 TSD 非空比例。
仅统计 Alu/L1/SVA。解读注意:
- polyA 面板是主判据,TSD 面板是"能确证的下限"——TSD 判定要求断点
  精确到碱基且两侧序列完全相同,真实 SV call 断点常有数 bp–数十 bp
  模糊,一侧偏移 1bp 即配对失败,检出率天然极低;
- Alu 的 polyA 检出率通常高于 L1:Alu 全长仅 ~300bp,插入序列常完整
  回收(含 3' 尾);L1 多为 5' 截短拷贝,caller 对长 ALT 的截断更激进;
- SVA 事件数通常极少,占比无统计意义时不下结论;
- "至少 X% 是 L1 介导",不是"只有 X% 是 L1 介导"(见判读不对称性)。

#### 四、参数与依据

| 参数 | 默认 | 依据 |
| --- | --- | --- |
| polyA 窗口 | 60 bp | 尾长典型 20–150bp;全序列统计会被长元件稀释 |
| 连续 A 阈值 | 12 | (1/4)^12≈6e-8 随机背景概率,特异性充分 |
| A 占比阈值 | 0.90 | 容忍测序错误打断的真实尾巴 |
| TSD k 范围 | 7–20 | TPRT 典型 TSD 长度分布(L1 ~15-20bp, Alu ~7-15bp) |
| TSD 侧翼 | 30 bp | 覆盖断点模糊度;匹配要求不随 flank 放宽 |
| polyA 最小插入 | 200 bp | 更短的插入多为测序噪声 |

#### 五、参考

- Feng Q, et al. Human L1 retrotransposon encodes a conserved endonuclease
  required for retrotransposition. *Cell* 87:905-916 (1996). — L1 ORF2p
  内切酶与 TPRT
- Luan DD, et al. Reverse transcription of R2Bm RNA is primed by a nick at
  the chromosomal target site. *Cell* 72:595-605 (1993). — TPRT 原始提出
- Moran JV, et al. High frequency retrotransposition in cultured mammalian
  cells. *Cell* 87:917-927 (1996). — L1 cis/preference 与反式动员
- Dewannieux M, Esnault C, Heidmann T. LINE-mediated retrotransposition of
  marked Alu sequences. *Nat Genet* 35:41-48 (2003). PMID 12897783 —
  Alu 被 L1 蛋白反式动员、3' poly(A) 对其移动必需的直接实验证据
  ("parasite's parasite" 之说出自同期评论 Nat Genet 35:15-16)
- Hancks DC, Kazazian HH. Active human retrotransposons: variation and
  disease. *Curr Opin Genet Dev* 22:191-203 (2012). PMID 22406018 —
  SVA/L1 插入特征综述(TSD 长度、polyA、5' 截短)

### 3.3 TRA / TE-TE 重组 —— 方法文档

日期: 2026-11-05 原独立文档并入本 README; 2026-11-23 拆出年轻插入验证部分(现 §3.4)
前置: §3.1(family 归桶规则)
本节只记录**方法**: 输入契约、算法步骤、判定规则、统计设计、参数依据、
质控点、方法学教训、局限。不含具体数据结果与运行用法; 结果见输出文件,
运行方式见代码 argparse 帮助与 sv.smk 的 sv_tra_te 规则。

---

#### 0. 问题定义与证据框架

**问题**:
- Q2 (TRA 线): 跨染色体易位是否由 TE 驱动(NAHR)? 依据什么分级判定?

**TRA 线证据漏斗(四级, 逐级收紧)**:

| 证据级 | 判定 | TE 驱动? |
|---|---|---|
| T0 事件 | 跨染色体 TRA(镜像去重后) | 未判定 |
| T1 相关 | ≥1 断点落 TE 内 | 否(单侧巧合可能) |
| T2 配对 | 双断点均落 TE 内, 且富集检验显著 | 候选 |
| T3 驱动 | 双拷贝序列同源达 NAHR 机制要求 | **是(序列级)** |

判定口径: 以跨染色体断点(inter-chromosomal TRA)为转座(易位)事件入口;
断点落 TE 注释为"TE 相关"; 两拷贝序列同源达到 NAHR 机制前提才判
"TE 驱动"。单侧 TE、双侧无同源、仅富集未过同源的事件, 各级分列,
不越级归因。TE 驱动的确定性: 双拷贝高同源是 NAHR 的机制前提本身
(重组必须配对同源序列), 统计富集 + 序列同源构成完整证据链;
直接实锤为接合序列本地组装(未系统化, 局限之一)。

---

#### 1. 输入契约

| 输入 | 关键字段/约定 | 用途 |
|---|---|---|
| SURVIVOR 合并 VCF | INFO: SVTYPE, CHR2, END, SUPP_VEC, STRANDS, CIPOS/CIEND | TRA/INS 事件来源 |
| 单样本 pbsv VCF | FORMAT=GT:AD:DP; BND 记录含 MATEID、ALT 括号方向 | VAF 验证、链方向恢复 |
| rmsk TE BED | chrom/start/end/name | 断点 TE 注释 |
| 参考 FASTA(.fai) | — | TE 拷贝序列提取 |
| 基因注释 GTF | gene 层, gene_name/gene_id | 断点基因注释 |
| fiber BAM | HiFi | 接合读长抽查 |

##### 1.1 坐标与 ID 约定(易错点)

- VCF POS: **1-based**; INS 插入位于 POS 与 POS+1 之间; BND 的 POS 是
  锚定碱基。
- BED: **0-based half-open**。断点判定: `bisect_right(starts, pos)-1`
  取候选区间, 验 `ends[i] > pos`(TE 长度 >> 1bp 边界误差, 1-based/
  0-based 边界差对判定无实质影响)。
- samtools faidx: **1-based inclusive**; BED→faidx 坐标转换: start+1。
- SUPP_VEC: 位数=样本数, 对照组在前。读取时 '01' 丢前导零变 '1',
  **必须补零对齐**(已踩过的 bug, 导致一组统计消失)。

##### 1.2 工具环境要点

- blastn: BLAST+ 参数是 `-num_threads`, **不是 `-threads`**(传错则
  进程报错退出、stdout 空);
- 容器内临时文件必须写绝对路径(apptainer 只 bind 数据盘, 容器内相对
  路径不可见);
- 外部工具调用必须检查 returncode 与输出非空。

---

#### 2. TRA 记录解析与独立事件判定

##### 2.1 解析

- 合并 VCF 取 SVTYPE=TRA; INFO 正则 `(\w+)=([^;]+)` 抽 CHR2/END(第二
  断点)/STRANDS/SUPP_VEC;
- 术语: pbsv 单样本报 BND, SURVIVOR 合并转为 TRA——同一事件在管线
  不同阶段的名字。TRA 标签 = inter-chromosomal BND 的合并态;
  intra-BND 属染色体内重排, 不在转座(易位)分析范围。

##### 2.2 镜像去重(记录级→事件级)

- 同一接合从两端各报一条(A→B 与 B→A), 去重 key =
  `sorted([(chrom1,pos1),(chrom2,pos2)])`, drop_duplicates;
- 质控: 记录数必须 = 2 × 事件数(每接合两端各一条, 无残缺);
- 局限: 仅合并坐标完全相同的镜像; 近邻断点仍各自计数(见 4.4)。

##### 2.3 链信息(SURVIVOR 有损压缩 —— 已实验定性)

迷你 VCF 实验(4 条覆盖全部括号方向 N[X[ / N]X] / ]X]N / [X[N 的 BND
记录跑 SURVIVOR merge): 输出**全部**被改写为 ALT=N[chr2:pos[ +
STRANDS=++; 遍历合并参数组合(type=0/1 × strand=0/1)均如此。

**结论**: SURVIVOR merge **没有任何参数可以保住链信息**。第 5 参数
(strand 0/1)只控制合并时是否按链分组, 不控制输出保真度; TRA 的
ALT/STRANDS 在写出时被硬编码重写, 方向信息确定性丢失。
**可恢复**: 单样本 pbsv VCF 保留原始括号方向(4 种模式齐全), merged
TRA 可按坐标回连单样本 BND, 真实链方向回查原始 ALT 即可恢复;
平衡/不平衡判定需回单样本 VCF + BAM 读深。

##### 2.4 复发断点簇(独立性评估)

- 全部断点按 (chrom, pos) 排序, 同染色体相邻间距 >2000bp 处切簇
  (cumsum(diff>2000));
- 簇内 ≥2 断点为复发簇;
- 判读: 簇横跨两组 → 群体复发热点或比对伪影, 非组特异获得性重排;
  近邻簇未合并时, 事件数可能高估。

---

#### 3. 断点 TE 富集检验(统计设计)

##### 3.1 断点-TE 判定

BED 按 (chrom, start) 排序建 per-chrom 数组, 二分查找判定(见 1.1)。
BED 列有混合类型脏数据, 数值化+去空后使用。

##### 3.2 三重对照设计(只报倍数没有说服力, 必须有零分布)

1. **基因组背景**: TE 总碱基 / 基因组总碱基(染色体长度用 max(end)
   近似, T2T 无 gap 时误差 <1%);
2. **随机置换-单断点**: N(默认 1000)次随机断点(染色体按长度加权,
   位置均匀), TE 命中率 → 用于 sanity check(应≈背景, 证明采样无偏);
3. **随机置换-双断点独立期望**: M(默认 2000)对独立随机断点对, 双
   命中率 → 零分布(≈单断点率², 符合独立假设);
4. **二项检验**: H0: 双 TE 率 = 独立期望; H1: >; 单侧 greater。
   置换与二项两种设计独立给出同一结论, 互为交叉验证。

##### 3.3 参数表

| 参数 | 默认 | 依据 |
|---|---|---|
| 置换次数(单断点) | 1000 | SE≈1.4%, 标准经验值 |
| 置换次数(双断点对) | 2000 | 同上 |
| 随机种子 | 42 | 可复现 |
| 显著性 | α=0.05 | 常规 |

---

#### 4. TE 家族配对与 shuffle 期望

##### 4.1 家族归类(前缀规则)

Alu / L1 / SVA / LTR/ERV / LINE2 / MIR / DNA / satellite-simple / other,
按 repeat 名前缀映射。**审核重点**: 手工映射表, 误归类直接影响期望值。
注意点: MacERV* 归 LTR/ERV; Alu 前缀(含 AluYRa1 等)归 Alu;
L1 开头(含 L1_RS2 等)归 L1; 名含 ERV 或 LTR/MLT/MER 前缀归 LTR/ERV。

##### 4.2 配对观测

both-TE 事件, 每事件两断点家族组成 pair = sorted(fam1, fam2) 拼接, 计数。

##### 4.3 shuffle 期望公式

设 both-TE 断点总数 N(= 2×事件数), 家族 f 频数 c_f。从 N 个断点随机
不放回抽两个为一对, 重复事件数次:

- 跨家族 (a≠b): E = 2·c_a·c_b / (N(N-1)) × 事件数
- 同家族 (a=b): E = c_a(c_a-1) / (N(N-1)) × 事件数

##### 4.4 判读规则

- ratio = obs/exp > 1.5 → 家族配对偏好(富集);
- ratio < 0.7 → 回避(低于随机);
- exp < 1 的小类别倍数不稳, 仅参考不作结论。

---

#### 5. 序列级同源验证(blastn)

##### 5.0 同源输入序列规格(序列背景取多大)

- **取区间**: 断点命中的**整条 rmsk TE 注释区间** [start, end)(0-based
  half-open), 即整条 TE 拷贝全长(代码: tra_te_analysis.py stage3,
  fetch_seq(starts[i], ends[i])); **不是**断点 ±N bp 固定窗口, 也**不含**
  任何侧翼序列。
- **长度**: 完全由 rmsk 注释决定(Alu ~300bp, L1 可达 ~6kb), 代码无长度
  上下限, blastn 输入长度随 family 波动 —— 不同 family 的 hom_index
  (pident × hsp_len) 不可直接横向比较。
- **重叠命中**: 断点落在多条重叠 TE 注释时, te_hit 二分只取第一条命中
  区间(见 1.1 bisect 规则), 不枚举全部重叠拷贝。
- **侧翼**: 同源臂只算 TE 体内; 真实 NAHR 同源臂延伸进侧翼时 hom_index
  低估(见局限汇总 TE 侧翼同源未纳入 条目)。
- **审查点**: NAHR 阈值(5.2: 90% / 100bp)的尺度语义依附于上述取法;
  若改为断点 ±N bp 或 min(全长, cap), 阈值依据须同步重估。

##### 5.1 流程(逐对执行)

1. rmsk BED 取两断点各自命中区间;
2. samtools faidx 提取参考序列(坐标转换见 1.1; 临时文件写绝对路径);
3. blastn 双链比对(outfmt 6, pident/length/坐标列);
4. 取最优 HSP(High-scoring Segment Pair,高分配对片段: blastn 输出的
   每个局部比对块即一个 HSP,含 pident/length/坐标;图中 "best HSP
   identity/length"、"no HSP" 均指它): max(pident × length)(匹配碱基
   数最多的高同源块);
5. 记录 ident(%)、hsp_len(bp)、hom_index = pident × hsp_len / 100。

##### 5.2 NAHR 判定阈值与依据

- **NAHR 级**: ident ≥ 90% 且 hsp_len ≥ 100bp。依据: 哺乳动物非等位
  同源重组效率要求 ~95%+ 同源、数百 bp 以上同源臂(经典 L1 近期拷贝对
  值); 90%/100bp 为保守下限, 留测序/比对噪声余量;
- **无同源**: blastn 无任何 HSP(默认 e-value)→ 排除 NAHR, 归
  NHEJ/巧合落点(断点碰巧落 TE——基因组 TE 覆盖高, 随机断裂常落 TE);
- 中间层(80-90%)单列, 不计入 NAHR。

##### 5.3 方法学教训(两次错误迭代, 引以为戒)

- v1 difflib SequenceMatcher 一致块: 只匹配正向相同序列, **不搜反向
  互补**——负链 TE 拷贝被误判不同源(错误结论: 同名对多数"不同源");
- v2 blastn 传参 `-threads`(正确为 `-num_threads`), 进程报错退出、
  stdout 空, 全部记 0(错误结论: "全部无同源");
- v3 修正后结果与生物学预期一致。
- 教训: (a) 双序列同源必须用真正比对器且搜双链; (b) 外部调用必须检查
  returncode 与输出非空; (c) **"全部为 0"是系统性错误信号, 不是生物
  学结论**; (d) 已知阳性对照必须在流水线中复现(质控点)。

---

#### 6. TRA 全景刻画方法

##### 6.1 染色体对网络

cpair = sorted(chrom1, chrom2) 计数; 热点判定按计数排序。

##### 6.2 基因注释

- GTF gene 层, per-chrom 按 gene start 排序建 bisect 数组;
- 断点判定: 候选 gene 附近 ±2 行内逐一验 `gene_start <= pos <=
  gene_end`(处理重叠基因); 记 gene_name(无则 gene_id);
- "断点在基因内" = 落 gene body(**不含**启动子窗口, 0 容忍);
- **融合候选** = 双断点均在 gene body 内。

##### 6.3 融合候选判读规则

仅断点层面注释; 无 RNA 数据时**不做转录层结论**; 小样本(n<10)下
组间计数失衡不作方向性结论(二项检验确认不显著后仅列候选)。

---

#### 7. 统计设计汇总

| 检验 | 位置 | H0 | 统计量 |
|---|---|---|---|
| 置换(单断点) | Stage1 | 断点随机 | N 次 TE 命中率(≈背景的 sanity check) |
| 置换(双断点) | Stage1 | 两断点独立 | M 对双命中率(零分布) |
| 二项(greater) | Stage1 | 双 TE 率 = 独立期望 | k/n 单侧检验 |
| shuffle 期望 | Stage2 | 配对随机 | 家族频数乘积期望 |
| 二项(小样本事后) | 全景 | 两组计数 1:1 | 失衡显著性 |

---

#### 8. 质控点清单(复现时逐项核对)

1. 置换单断点率 ≈ 基因组背景(采样无偏 sanity check);
2. blastn 已知阳性对照(手工验证过的高同源对)在流水线 stage3 复现;
3. 镜像去重计数: 记录数 = 2 × 事件数;
4. SUPP_VEC 补零后三组计数与 VCF 直接 grep 一致。

---

#### 9. 局限汇总(按影响排序)

1. **1v1 设计**: "组特异"= 个体私有, 无年龄效应统计力; 一切方向性
   结论仅候选级;
2. **接合序列未组装**: NAHR 判定靠拷贝同源度 + 读长比对签名, 非接合
   序列直接证据; 严格实锤需逐事件本地组装;
3. **STRANDS 信息丢失(可恢复)**: SURVIVOR 写出时确定性丢链信息(无
   参数可解, 见 4.3 实测); 可按坐标回查单样本 VCF 恢复;
4. **独立事件数可能高估**: 近邻断点簇未合并(可能是同一复杂重排的
   多次拆分);
5. 融合候选无转录层验证;
6. TE 侧翼同源未纳入 NAHR 同源臂计算(hom_index 低估)。


### 3.4 年轻插入验证 —— 方法文档

日期: 2026-11-23 自原 TRA / TE-TE 重组与年轻插入验证 文档拆出(MEI 线部分)
前置: §3.1(MEI 管线总览) / §3.2(polyA+TSD 原理)
问题: Q1 (MEI 线): INS 中哪些是转座事件? 是否有新发(体细胞)转座?
本节给年轻/新发判据: 转座判定(INS 构成分类)+ 定年筛选 + VAF 验证;
TRA 证据漏斗与序列级同源验证见 §3.3。运行入口: 定年
bin/utils/te_age_analysis.py; 事件表/定年表由 bin/mei_analysis.py 产出。

---

#### 0. 输入契约与 id 映射

| 输入 | 关键字段/约定 | 用途 |
|---|---|---|
| MEI 注释表 (mei_events.tsv) | family/polyA/tsd/context/gene/te_repeat/in_te_region | 事件分类输入 |
| 定年表 (te_divergence_per_ins.tsv) | w_div(Kimura), jc_div, mu 档年龄 | 候选筛选 |
| merged VCF / 参考 FASTA | 见 §3.1 输入契约 | INS 来源、TSD 断点侧翼 |

##### 0.1 ins.fa 与事件表 id 映射(易错点)

- MEI 表与 ins.fa 的 id 映射: ins.fa 第 i 条序列 = 事件表中超过长度
  阈值的第 i 条(VCF 记录顺序); 需校验 best-hit repeat 名一致(0 错配
  才能继续)。
---

#### 1. INS 构成分类(转座 vs 非转座)

**方法**: 事件表 family 列计数。family 由 RepeatMasker 分类产出:
插入 ALT 序列超阈值者写 FASTA → RepeatMasker `-s -no_is -lib <te_lib>`
→ 每条取**覆盖跨度最大**的命中(非 SW 最高)→ 按 repeat 名前缀归桶:

| family 桶 | 规则 |
|---|---|
| Alu | 名以 Alu 开头 (SINE) |
| L1 | 名以 L1 开头 (LINE) |
| SVA | 名以 SVA 开头 |
| LTR | 名含 ERV 或以 LTR/MLT/MER/THE 开头 |
| DNA | 名以 Tigger/hAT/Charlie/Zaphod 等前缀开头 |
| simple_repeat / microsatellite | 名含 (X)n 式括号模式 |
| short_unclassified | 长度低于阈值且无命中 |
| other / unclassified | 卫星等其余; 未命中且超阈值 |

**判定逻辑三层区分**(详见 MEI_METHODS.md 第三节): 序列来源层
(RepeatMasker 同源)≠ 移动机制层(TPRT / 整合酶 / 剪切-粘贴)≠
事件层(新生插入 / 历史残留 / 伪影, 需指纹+定年证据分级)。
simple/microsatellite = 串连重复扩张, 非转座。
未分类处理**声明为假设**(假设非 TE → 下界; 全 TE → 上界),
结论以区间表达, 不给单点值。

---

#### 2. VAF 验证方法(生殖系多态 vs 体细胞嵌合)

##### 2.1 候选筛选(确定性规则)

- 事件表与定年表按 (chrom, pos) 等值 merge(同源 VCF 坐标, 无偏移);
- 候选 = (family ∈ {Alu,L1,SVA} 且 w_div < 1.0) ∪ (context='exon' 且
  w_div < 2.5);
- (chrom, pos) 去重。

##### 2.2 原始支持数回查

合并 VCF 丢失 FORMAT(只剩 SUPP_VEC), 必须回查单样本 pbsv VCF:

1. 全量扫描建索引 `CHROM_POS -> FORMAT 字段原文`;
2. 候选 POS 直接查表(合并 POS 与单样本 POS 同源; 若发现偏移改用
   ±1bp 容差);
3. 解析 GT:AD:DP, AD = (ref, alt)。

##### 2.3 判定规则

- VAF = alt / (ref + alt);
- 生殖系杂合: GT=0/1 且 VAF≈0.5(抽样噪声下不设硬阈值, 报中位数与
  范围); 纯合: GT=1/1 且 VAF=1.0;
- **体细胞嵌合判据**: VAF 显著低于 0.5(经验 <0.25)且深度足够;
- 质控: 回查的组别归属须与 SUPP_VEC 自洽(三方交叉验证 merge 正确性);
  深度需足以区分 0.5 与低 VAF。

---

#### 3. 统计设计(定年分布比较)

| 检验 | 位置 | H0 | 统计量 |
|---|---|---|---|
| MWU/KS | te_age | 组间年龄分布相同 | div 分布差异 |

---

#### 4. 质控点

1. 事件表 ↔ ins.fa 一致性: 长度过滤后条数相等, best-hit repeat 名
   0 错配;
2. VAF 回查组别 vs SUPP_VEC 自洽(三方一致);

---

#### 5. 局限

1. VAF 取自 caller 的 AD, 未回 BAM 独立复核跨断点读数;


## 4. 参考文献

1. 方法文档: 本 README §3(与代码同版本, 含参数依据与方法学教训)。
2. Ensembl VEP —— https://www.ensembl.org/info/docs/tools/vep/index.html
3. RepeatMasker —— https://www.repeatmasker.org;Dfam —— https://www.dfam.org
4. COSMIC —— https://cancer.sanger.ac.uk/cosmic(oncoprint 基因列表来源)
5. 运行环境:`sv.yaml` / `sv.def`(sv.sif:RepeatMasker/blastn/pysam/
   pandas/matplotlib);samtools 用独立 SIF
   (`/data/16T/Database/env/samtools/samtools.sif`)。
