# MEI 指纹检测原理:poly(A) 尾与 TSD —— 为什么它们是 L1 介导插入的标记

配套脚本:`modules/sv/bin/mei_analysis.py`(检测) + `mei_plot.py` polya_tsd_fingerprint 图(可视化)。
管线总览见 `MEI_METHODS.md`;本文档聚焦两个指纹的生物学原理与算法细节。

---

## 一、背景:MEI 与 L1 介导的逆转座

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
mei_analysis.py 只对 family ∈ {L1, Alu, SVA} 做 TSD 检测、polya_tsd 图只画这三个
家族的原因。

## 二、TPRT:指纹的产生机制

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

## 三、算法实现(modules/sv/bin/mei_analysis.py)

### 3.1 poly(A) 尾检测 — `poly_a_tail()`

```python
def poly_a_tail(seq, min_len=12, min_frac=0.9):
    tail_win = seq[-60:] if len(seq) >= 60 else seq
    best_run = max连续A数(tail_win)          # 允许被打断,取最长纯A run
    a_frac   = count(A)/len(tail_win)
    return best_run >= min_len or a_frac >= min_frac
```

- 窗口:插入序列**末 60 bp**(不足 60 取全长)。poly(A) 尾长度典型
  20–150 bp,60 bp 窗口兼顾灵敏度与特异性(全序列统计会被长元件稀释)。
- 判阳条件(满足其一):
  1. 窗口内**连续 A ≥12 bp**——真实测序/比对错误率下,随机序列出现
     12 连 A 的概率 ≈ (1/4)^12 ≈ 6e-8,基本排除背景;
  2. 窗口内 **A 占比 ≥90%**——容忍测序错误造成的零星打断,覆盖
     "A-rich 但不纯"的真实尾巴(逆转录/复制过程中的滑移、替换)。
- 前置条件:`ins_len ≥ 200` 才检测(太短的插入多为噪声,记 NaN)。
- 只查 3' 端、不区分链方向:MEI 的 ALT 序列方向由 caller 决定,polyA
  若在序列 5' 端说明该序列是反向的,靠 RepeatMasker 的 strand 列交叉
  判断(见 mei_events.tsv 的 te_strand)。

### 3.2 TSD 检测 — `microhomology_tsd()`

```python
def microhomology_tsd(ref_fa_idx, chrom, pos, ins_seq, flank=30):
    up   = ref[pos-1-30 : pos-1]     # 断点上游 30bp(VCF POS 为 1-based 锚定碱基)
    down = ref[pos   : pos+1+30]     # 断点下游 30bp(插入位于 POS 与 POS+1 之间)
    for k in range(20, 6, -1):       # 从长到短,优先报最长 TSD
        a, b = up[-k:], down[:k]     # 紧贴断点两侧的 k-mer
        if a == b and set(a) != {'N'}:
            return {'tsd': a, 'len': k}
    return None
```

- 坐标约定:VCF 中 INS 的 POS 是插入前 1-based 锚定碱基,插入发生在
  POS 与 POS+1 之间;pysam fetch 是 0-based 半开区间,故上游取
  `[pos-1-flank, pos)`,下游取 `[pos, pos+flank+1)`。
- 判定标准:紧邻断点上游末 k bp 与下游首 k bp**完全相同**(k=20 递减到
  7)。要求"紧邻 + 两侧同时存在 + 精确匹配",对应 TPRT 产生的典型
  7–20 bp 直接重复;全 N(参考缺口)排除。
- 从长到短扫描:一个真实 TSD 若两侧完好,最长匹配就是它本身;先扫短
  k 会把 TSD 的子串误当 TSD,长度信息也丢失。
- 只对 family ∈ {L1, Alu, SVA} 检测(LTR/DNA 转座子无此机制,检测无意义,
  且会引入假阳性)。
- 输出:mei_events.tsv 的 `tsd` 列存重复序列本身,`tsd_len` 存长度,
  未检出为空。

### 3.3 下游统计(mei_plot.py fig6)

- polyA 率 = family 内 `polyA == True` 的比例;
- TSD 率 = family 内 `tsd` 非空的比例;
- 仅统计 Alu/L1/SVA,以"携带指纹的插入占比"作为该批事件中 L1(TPRT)
  介导新生插入比例的**下限估计**。

## 四、如何解读 fig6(以 luanchao_mulatta youth_vs_aged 数据集为例)

实测值(mei_events.tsv 复核):Alu n=45,543,L1 n=5,587,SVA n=20。

| 家族 | polyA 率 | TSD 率 |
| --- | --- | --- |
| Alu | 40.4% | 1.2% |
| L1  | 21.4% | 0.9% |
| SVA | 0.0% | 0.0%(n=20,无统计意义) |

1. **polyA 面板是主要证据**:约 1/5–2/5 的 Alu/L1 插入 3' 端带 A-rich
   尾,直接支持这部分事件是 TPRT 介导的新生插入,而非比对错误或参考中
   已有拷贝的误报。
2. **Alu 检出率高于 L1 的原因**:Alu 全长仅 ~300 bp,插入序列通常完整
   回收(含 3' 尾);L1 插入多为 5' 截短拷贝,caller 对长 ALT 的表示/
   截断也更激进,尾巴更易丢失。
3. **TSD 检出率极低是分辨率问题,不是机制缺失**:TSD 判定要求断点精确
   到碱基且两侧序列完全相同;pbsv/SURVIVOR 的 INS 断点常有数 bp–数十 bp
   模糊,一侧偏移 1 bp 即配对失败。本数据集检出的 TSD 几乎全为
   poly(A)-like 同聚物(如 AAAAAAAAAA),也说明精确匹配在真实 SV call
   上何等苛刻。**TSD 面板应视为"能确证的下限",polyA 面板才是主判据。**
4. **未检出 polyA 的 L1/Alu 事件不能就此否定**:ALT 截断(见
   MEI_METHODS.md 已知限制 #2)会人为压低检出率,故 fig6 的读法是
   "至少 X% 是 L1 介导",而不是"只有 X% 是 L1 介导"。

## 五、参数与调优

| 参数 | 默认 | 位置 | 调整建议 |
| --- | --- | --- | --- |
| polyA 窗口 | 60 bp | `poly_a_tail()` 硬编码 | 想抓更长的尾巴可加大,但特异性略降 |
| 连续 A 阈值 | 12 | 同上 | 降 to 10 会略增灵敏度,(1/4)^10≈1e-6 仍可控 |
| A 占比阈值 | 0.90 | 同上 | 测序质量差的数据集可降到 0.85 |
| TSD k 范围 | 7–20 | `microhomology_tsd()` 循环 | 物种/家族 TSD 长度分布不同时可收窄 |
| TSD 侧翼 | 30 bp | 同上 flank | 断点误差大的 caller 可加大,但匹配要求不变 |
| polyA 最小插入 | 200 bp | 主循环 | — |

## 六、参考

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
