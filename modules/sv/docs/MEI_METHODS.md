# MEI (Mobile Element Insertion) 分析方法文档

针对 SURVIVOR 合并 SV VCF(pbsv 来源)的转座事件分析管线的方法说明。
只记录方法(输入契约、算法、判定规则、参数依据、局限), 不含具体数据结果
与运行用法; 运行方式见代码入口的 argparse 帮助与工作流调用
(`modules/sv/sv.smk` sv_mei 规则)。
代码入口: `modules/sv/bin/mei_analysis.py`(唯一入口; 绘图经
`bin/utils/mei_plot.py` 串联)。TE 插入定年(divergence-from-consensus
分子钟)为独立分析(`bin/utils/te_age_analysis.py`), 不在本管线内。
poly(A)/TSD 指纹的检测算法与 TPRT 生物学原理详见 `README_polyA_TSD.md`。

## 一、输入契约

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

## 二、分析流程(五步)

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

## 三、判定逻辑(三层区分, 勿混用)

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

## 四、输出 schema(接口定义, 非用法)

per-event 表(坐标/长度/组别/家族/亚家族/覆盖度/polyA/TSD/语境/基因/
TE 区域/全长 L1 标记)、家族×组别计数矩阵、语境计数、全长 L1 清单、
汇总 json、图组(family_by_group / insert_length_by_family /
genomic_context / te_region_ratio / full_length_L1 /
polya_tsd_fingerprint)、自包含 HTML 报告、RepeatMasker 原始输出留档
(可复查亚家族)。具体文件名见代码 docstring。

## 五、方法学局限

1. **1v1 设计**: 每组单样本时"组特异"是个体多态而非组效应; 组间定量
   比较需每组多重复(SUPP_VEC 逐位自适应, 方法不变)。
2. **polyA 检出率是下限**: caller 对 ALT 序列有截断, 3' 端常缺失,
   人为压低检出率; 长元件(如 L1, 多为 5' 截短拷贝)受影响更大。
3. **TSD 仅阳性可确证**: 参考位点已有祖先插入时断点两侧无重复序列;
   且断点 1bp 偏移即配对失败, 真实 SV call 上检出率天然极低——
   TSD 面板是"能确证的下限", polyA 才是主判据。
4. **RepeatMasker 分类局限**: 比对的是共识序列而非基因组拷贝; 嵌合/
   重排插入只报最长命中; 低覆盖度命中分类可信度低。

## 六、差异分析(mei_diff, 探索性)

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
