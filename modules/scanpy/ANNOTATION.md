# scRNA-seq 自动注释与双向证据体系

本文档描述 `modules/scanpy` 中 `scRNAseq.py` 的自动注释逻辑(`mode_auto`)以及配套的
双向证据(two-sided evidence)审核体系。内容对应当前代码实现,修改代码时同步更新本文档。

数据流位置:`scanpy_merge(merged.h5ad)` → **`scanpy_auto(auto.h5ad, 注释)`** →
`scanpy_reannotate(人工覆写,可选)` → `scanpy_advanced` / `scanpy_pseudobulk`。

---

## 1. mode_auto 自动注释逻辑

### 1.1 Step 0:组织面板(tissue panel)

注释前先向 LLM 查询组织的已知细胞类型(`_query_tissue_cell_types`):

```
tissue_cell_types = {
  "Stromal fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL15A1", "COL6A1"],
  "Myofibroblast": ["ACTA2", "TAGLN", "MYL9", "CALD1", "COL1A1", "COL3A1", "PDGFRB", "DES"],
  ...
}   # Uterus 面板为 34 类
```

面板有三重作用:

1. **命名天花板**:提示词强制使用面板中的精确名称(允许在 marker 有诊断性时使用文献亚型名)。
   面板粒度直接决定注释粒度——面板里基质只有一类时,7 个基质亚状态会被压成一个标签。
2. **canonical marker 来源**:每类的标准 marker 是"身份签名",区别于 DEG 证据。
3. **双向证据原料**:竞争类型的 marker 集用于差分/排除论证(见 §3)。

面板持久化在 `audit_report.json` 的 `tissue_cell_types` 字段。

### 1.2 迭代循环

```
while iteration < max_iterations:
    Step 1   _auto_cluster_iteration   normalize → HVG → PCA → harmony → leiden → DEG
             (迭代>1 时先从 raw_adata 恢复全基因原始计数)
    Step 2   _auto_annotate_all_clusters  逐 cluster LLM 注释
    Step 2.1 _auto_normalize_annotation_names  标签对齐面板命名
    Step 2.5 _auto_refine_annotations  特异性加权精修
    Step 3   _auto_quality_reports      质控判定(应过滤的 cluster)
    Step 3.5 _check_cluster_continuity  空间连续性(断裂→升分辨率重聚)
    Step 3.75 _check_cell_type_separation  同型分离度(重聚建议)
    Step 4   _filter_flagged_cells      过滤 flagged 细胞 → 重聚;或干净→退出
```

`raw_adata`(全基因原始计数快照)在整个循环中不可变,是每次重聚的起点;
最终输出 h5ad 的 `layers["counts"]` 为全基因整数计数(pseudobulk 的计数来源)。

### 1.3 单 cluster 注释上下文(prompt 输入)

按证据权重排序:

1. **PRIMARY:特异性表**(`_compute_marker_specificity`,面板 marker × 本 cluster):
   `specificity = this_mean / other_max`,附 pct 表达。≥2.0x 为强证据。
2. **PANEL MARKER EVIDENCE - COMPLETE(双向证据表)**:面板全部 canonical marker 的
   `gene (mean, %, spec)` 三元组,按细胞类型分组。这是 presence/absence 断言的权威依据
   (见 §3.5 的失败教训)。
3. **SECONDARY:Top 50 DEGs**(Wilcoxon 排名)——统计显著性≠细胞类型特异性,只作补充。
4. 其他 cluster 的已有注释 + UMAP 距离、既往错注反馈(misannotation history)、组织面板。

**断言契约(prompt 硬约束)**:

- `% 高 + spec 低` → "present, not specific",禁止写 "absent";
- `% <5%` → 可断言缺失,且 reasoning 必须引用该 % 值;
- 表中没有的基因 → 写 "not evaluated",禁止断言缺失。

### 1.4 LLM 输出 schema

```json
{
  "cell_type": "...",
  "key_markers": ["..."],          // 本 cluster 表达谱证据(DEG 子集),非身份签名
  "canonical_markers": ["..."],    // 文献标准 marker,身份签名
  "evidence": {                    // 组合型身份必填;单基因诊断型省略
    "programs": {"contractile": ["ACTA2","TAGLN","CALD1","PDGFRB"],
                 "ECM": ["COL1A1","COL3A1","COL4A2","COL18A1"]},
    "ruling_out": {"Smooth muscle cell": ["MYH11","DES"],
                   "Stromal fibroblast": ["DCN","LUM"]}
  },
  "reasoning": "...",
  "confidence": "high|medium|low",
  "quality_flag": null|"low_quality"|"unannotated"|...
}
```

### 1.5 后处理与产物

- `Uterus_scTE_auto.h5ad`:obs 含 `leiden / cell_type`;`uns["cluster_annotation_overrides"]`
  为覆写审计轨迹;`layers["counts"]` 全基因整数。
- `Uterus_scTE_auto_reports/table/`:`annotation_report.tsv`(每 cluster 一行,
  confidence/key_markers/flags/should_merge/reasoning)、`annotation_audit.json`
  (独立表达核验 + 面板 + 迭代上下文 ctx,是注释元数据的**权威副本**)、
  `audit_report.json`、`references.tsv`、`decision_log.md`。
- `plot/marker/<cell_type>.png`:canonical marker UMAP 面板(见 §3.4)。

---

## 2. 双向证据体系

### 2.1 原则

单基因特异性不存在时(基质/成纤维谱系常态),细胞身份由**双向证据**定义:

- **正向**:"我表达什么"——诊断性 marker(diagnostic)或程序组合(composite)。
- **反向**:"我不表达什么"——竞争身份的标志基因缺失(ruling out),阴性证据同样可定量。

图的形式必须与证据的形状对齐——不同证据类型用不同审核图,不再一刀切画 marker UMAP。

### 2.2 证据类型

| 类型 | 判定 | 含义 | 审核图 |
|---|---|---|---|
| **diagnostic** | top marker ≥2x 特异 | 单基因即可定位 | canonical marker UMAP 面板 |
| **composite** | 无 ≥2x marker | 身份=程序组合(AND) | 签名图包(§3.3) |
| **negative** | 排除法为主 | 靠竞争身份标志缺失 | 差分 dotplot(§3.4) |
| **gradient** | 兄弟簇高重合 | 梯度/亚状态存疑 | 兄弟簇对照,倾向收敛父级 |

粒度铁律(C13 案例确立):**诊断性 marker 不足就保持父级标签**;
反向纠偏(C6 案例)——标签与区室/签名矛盾时改名,与加细粒度是两回事。

### 2.3 marker 源分级(图上标题的后缀即出处)

1. `canonical markers` — 组织面板/Literature 标准签名(身份验证用,首选)
2. LLM 返回的 per-cluster `canonical_markers`
3. `key markers (DEG)` — 表达谱证据,仅兜底;画出来"到处都红"是正常的,它不是身份签名
4. `+ supplement` — 人工补充的判别 marker(`-X/--marker-supplement`,见下),与上者叠加显示

**marker supplement(判别 marker 补充)**:近缘细胞类型(如 myofibroblast /
smooth muscle / stromal fibroblast)的 canonical 面板高度重叠,单靠面板无法判读
特异性。补充 JSON 按细胞类型给出组合判别 marker:

```json
{
  "Myofibroblast": {"positive": ["POSTN", "CTHRC1"],
                    "negative": ["MYH11", "MYOCD", "PI16", "DPT"]}
}
```

`positive` 期望表达,`negative` 期望不表达(对竞争类型的排除证据),后者在面板
标题带 `(exp. absent)` 后缀。键与细胞类型名大小写不敏感匹配,矩阵中缺失的基因
记 warning 后跳过;`_` 开头的键为注释键,忽略。接入方式:
`Params.scanpy.<counter>.auto.marker_supplement = {tissue: json 路径}`(scanpy.smk
在 scanpy_auto / scanpy_reannotate 两个规则透传 `-X`)。**仅作用于绘图**
(marker panels),不进入 LLM 提示词与特异性打分。来源有两路,按 cell_type key
合并,人工 `-X` 条目整条覆盖 LLM 自动版(大小写不敏感):
1. **LLM 自动产出**(默认开,`marker_supplement_llm: false` 或 CLI `-M` 关闭):
   auto 模式注释收敛后追加一次 LLM 调用,为近缘易混类型生成判别组合
   (positive 期望表达 + negative 期望不表达),原子写
   `<output stem>_reports/table/marker_supplement.auto.json`;reannotate
   自动发现同目录文件并合并。任何失败(调用/解析)降级为空,不阻塞运行。
2. **人工维护**:`Params.scanpy.<counter>.auto.marker_supplement` 的 JSON
   (需要覆盖或补充 LLM 版时写)。
示例:`modules/scanpy/examples/marker_supplement.example.json`
(真实案例:`/data/16T/luosg/data/Uterus_scTE_marker_supplement.json`)。刷新面板用 no-op 覆写
(JSON 中 corrected == original)+ `--forcerun scRNAseq_scanpy_reannotate
--rerun-triggers mtime`:figure-only rerun,只重建图与报告,h5ad 不动(§3.2)。

### 2.4 审核图包

**diagnostic → `plot/marker/<ct>.png`**
UMAP 特征图面板(canonical 集)+ Cluster ref 尾图(蓝色标目标簇)。
多 cluster 共享标签时合并高亮,标题 `[C0, C6, ...] — canonical markers`。

**composite → `plot/evidence/<ct>.png`(签名图包)+ `table/<ct>_signature.tsv`**

面板序列:

1. `negative control (random n=N)` — 同流水线跑随机基因集(seed 42):校准地板,
   排除"高表达区自然发亮"伪影;
2. `signature: min(z) AND-score` — **组合分数 = min(z_程序1, z_程序2, ...)**:
   AND 逻辑,单程序高不加分。禁止用 union 均值(那是 OR 逻辑,单个高表达程序会拉满,
   C6 案例教训);
3. `program z: <name>` — 逐程序 z 分数,展示各程序单独的分布(解释为什么单程序不特异);

各面板均带目标簇红圈凸包 + 编号标注;独立 `Cluster ref` 尾图(与 marker 图一致的簇边界
对照)待补,见 §5。

图注写明公式:`score = scanpy score_genes(表达分箱对照基因) → 全细胞 z 标准化;
combined = min(z) AND 逻辑;control = 同管线随机基因`。

`<ct>_signature.tsv`:每 cluster 的 pct 双阳性(排除混合物嫌疑)、各程序/combined/control
的 mean z——**图读数为辅,数字为准**,combined 的 argmax cluster 可从表中直接核对。

**negative → `plot/evidence/<ct>_differential.png` + `table/<ct>_evidence.tsv`**
差分 dotplot:行=排除 marker(按竞争身份分组),列=全部 cluster;
颜色=按基因行内标准化均值(行内相对比较),点大小=表达比例。
阴性行的读法:目标簇暗、竞争簇亮。TSV 为每基因×每 cluster 的 mean/pct 精确值。

### 2.5 双向证据在提示词中的实现

- 正向:特异性表 + 完整面板证据表(§1.3);
- 反向:`evidence.ruling_out` 声明排除性 marker;断言契约保证阴性断言可核验;
- 机械判定与 LLM 自报解耦:证据类型由管线按 specificity 机械判定,LLM 的
  `evidence` 块只提供结构(programs/ruling_out)。

---

## 3. 人工覆写与 scanpy_reannotate

### 3.1 覆写 JSON 契约

路径:`Params.scanpy.<counter>.auto.cluster_annotations.<tissue>`(per-tissue 映射):

```json
{
  "<cluster_id>": {
    "original_cell_type": "Myometrial fibroblast",
    "corrected_cell_type": "Myofibroblast",
    "canonical_markers": ["..."],            // 可选,标签不在面板时供图
    "evidence_markers": {"negative": [...], "positive": [...]},  // 可选,差分图
    "evidence": {"programs": {...}, "ruling_out": {...}},        // 可选,签名图包
    "reasoning": "..."
  }
}
```

完整示例(两种条目:纯改名、带 evidence 全字段):
`modules/scanpy/examples/cluster_annotations.example.json`。注意与 marker_supplement
不同,顶层键必须是 cluster id,**不支持注释键**(`_comment` 会触发 unknown cluster
id 报错)。

### 3.2 安全闸与就地覆盖语义

- **original_cell_type 必须匹配当前标签**(幂等重放例外:已是 corrected 则放行),
  不符即 ValueError——防止重聚类后 leiden 编号漂移导致误改(fail loudly)。
- **就地覆盖**:覆写后 `<auto>.bak` 保存原始版(**只写一次**,后续重跑不覆盖),
  原子替换(auto.h5ad → tmp → os.replace),报告/图在 auto 自己的目录内重生成,
  无平行 `_reannotated` 产物集。
- **幂等/增量**:标签无实际变化时跳过 h5ad 重写与审计重生成(纯图迭代只付画图成本);
  报告重建从 `audit_report.json` 的 ctx 继承原始 reasoning(报告只能从权威源派生)。
- categorical 收紧:覆写后 `remove_unused_categories` + 调色板按名对齐
  (否则 ghost 图例/颜色错位)。
- marker 图按 §2.3 源分级;改名 cluster 可通过 override 字段带签名/证据图。

---

## 4. 实战案例(设计的来源)

| 案例 | 问题 | 确立的设计 |
|---|---|---|
| **C13 Venous→Endothelial** | 2 个 <2x 静脉 marker 支撑细注 | 粒度铁律(§2.2);兄弟簇对照 |
| **C6 Stromal→Myofibroblast** | 无任何 ≥2x marker;LLM 断言 "absence of ACTA2" 与数据(82%)矛盾 | 双向证据体系整体;断言契约;AND 分数;面板差分 |
| key_markers 图泛红 | DEG 证据被当身份签名画 | marker 源分级(§2.3) |
| 覆写后图例残留 Venous | categorical 残留类目 | 覆写后类别清理(§3.2) |
| reasoning 占位语覆盖原文 | 报告重建丢失 LLM 证据链 | audit ctx 为权威源(§3.2) |

**LLM reasoning 的定位**:假设,不是证据。一切断言(尤其阴性断言)转述前必须对
表达矩阵核实;reasoning 只提供核查线索。

---

## 5. 已知缺口 / 待办

- [ ] 签名图包待补:**Cluster ref 尾图**(全簇边界对照,与 marker 图末图一致)。
- [ ] 签名图包待补:**表达深度参考面板**(每细胞 total UMI / n_genes)——排除
  "C6 整体高表达"导致的深度混淆(审阅者需对比 signature 与深度图是否同形);
  TSV 加每 cluster 深度列。
- [ ] `_compute_marker_specificity` 的 pct/均值来自 target cluster,但断言契约的
  机械核验(解析 reasoning 中的 absence 断言逐条对账)尚未实现——目前靠提示词约束。
- [ ] LLM 端到端(evidence 块 + 完整面板表)未在真实 auto 运行中验证
  (prompt 构建层已单元验证)。
- [ ] `should_merge` 标志(如 C21→C1、C18→C7)无人消费——待决定落实合并或输出
  suggested-merges 清单。
- [ ] 梯度型(gradient)审核图(兄弟簇对照)未实现,当前依赖人工判读。
