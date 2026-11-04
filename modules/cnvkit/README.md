# CNVkit 模块

基于 [CNVkit](https://cnvkit.readthedocs.io/) 的拷贝数变异(CNV)检测,从 gatk_prepare 的
markdup BAM 出发,产出每个样本的 bin 级 log2 比值(`.cnr`)与片段级分割(`.cns`)。

## 三种检测模式

CNVkit batch 的参考(reference)决定 log2 的参照系。本模块用
`Params.cnvkit.mode` 显式选择,`auto`(默认)由 node.py 从实验设计推导:

| mode | 参考 | 规则 | auto 推导条件 |
|---|---|---|---|
| `paired` | 每个 tumor 用**自己的 matched normal** 单独建参考 | `cnvkit_paired`(每对一个 job) | meta design 列存在配对(ctrl/exp) |
| `pooled` | 一批 normal 样本池化建共享参考(PON) | `cnvkit_reference` + `cnvkit_batch --reference` | 无配对但 `Params.workflow.control_samples` 非空 |
| `flat` | 无参考,假设期望覆盖均匀 | `cnvkit_batch` + 裸 `-n` | 既无配对也无 control_samples |

显式值优先于 `auto`;矛盾组合在 node.py 硬报错:
`paired` 无配对、`pooled` 无 control_samples 都会直接失败,
推导为 `paired` 时 `control_samples` 被清空(配对参考与池化参考不并存)。

## 各模式优劣(算法依据)

算法细节——biweight 池化、fix 偏差校正、可辨识性定理、单 normal 参考的噪声
行为——见源码级解析文档:
`CNVkit算法原理_源码级解析.md`
(以下 § 引用均指该文档)。

### paired(配对模式)

每个 tumor 独立执行 `cnvkit batch tumor.bam -n matched_normal.bam`,
用自己的 normal 建单样本参考。

**优势**:

1. **胚系 CNV 自动抵消**:肿瘤与 matched normal 携带同一套胚系缺失/扩增,
   log2 相减后只余体细胞事件——无需事后 DGV 过滤,这是配对模式最大的价值;
2. **批次效应抵消**:配对样本通常同批建库测序,文库制备/捕获的系统性偏差
   在相减中一同消掉;
3. **可辨识性成立**(§2.1):matched normal 对"体细胞 CNV"而言是 signal-free
   数据,扣除干净在数学上合法。

**劣势**:

1. **参考是 n=1,无平均**:pooled 参考的 biweight 中心(§2.2)把 N 个 normal
   的 bin 噪声平均掉,N 越大基线越稳;单 normal 参考把该样本自身的采样噪声、
   GC/捕获偏差的估计误差全部带入每个 bin 的期望值 → log2 分布更散、分割更易碎;
2. **跨样本不可比**:每个 tumor 对自己的 normal 归一,log2 之间没有共同标尺,
   队列级统一阈值/聚类要谨慎;
3. **算力重复**:多个 tumor 共享同一 normal 时(如 E14-1 配 4 个 tumor),
   该 normal 的覆盖被重算、参考被重建多次——结果正确,只是浪费。

**适用**:有 matched normal、目标是体细胞 CNV(肿瘤 vs 配对癌旁/血液)。

### pooled(基线/PON 模式)

`cnvkit_reference` 把 `control_samples` 的 BAM 一起喂给 batch,池化建
`reference.cnn`;所有 tumor 对这条共享基线调用(§2.2 的 biweight 池化)。

**优势**:

1. **噪声平均掉**:N 个 normal 的 bin 噪声按 1/√N 收缩,log2 在正常区域更紧地
   贴 0,分割干净;spread 列(§2.3)还让分歧大的 bin 自动降权(§3.4 fancy_wt);
2. **跨样本可比**:所有 tumor 归一到同一参照系,队列级比较、heatmap、统一阈值
   都成立;
3. **一次建参考反复用**:大队列只付一次建参考的成本。

**劣势**:

1. **胚系 CNV 混入报告**:panel 里没有患者的胚系变异,tumor 携带的多态性
   缺失/扩增会作为事件出现在结果里(胚系+体细胞混合清单),需事后过滤;
2. **批次敏感**:panel 与 tumor 不同批时,批次伪影变成假阳性(可用 batch 的
   `--cluster` 缓解,§2.3);
3. **panel 质量决定一切**(§2.6):PON 必须是匹配组织类型的正常样本,混入
   肿瘤样本等于往基线里掺信号,calling 向"混池平均图谱"收缩;官方经验
   panel ≥3 个 normal 才稳定,n<3 时可能不如 paired。

**适用**:无 matched normal 但有一批正常对照(cell line、cfDNA 队列、共享
panel 的项目);需要队列级一致标尺的场景。

### flat(无参考)

`cnvkit_batch` 传裸 `-n`,CNVkit 构建 flat reference——假设期望覆盖全基因组均匀。

**优势**:

1. **零外部依赖**:不需要任何 normal,单样本即可跑;
2. GC/重复序列协变量仍从 fasta 计算,fix 的滑窗偏差校正(§3.3)仍会执行。

**劣势**:

1. **校正最弱**:均匀期望与真实捕获偏差的差距全靠样本自身的滑窗中位数扣除,
   而 §2.6 的定量对照表明,样本自身回归/滑窗在校正同时会吃掉与协变量相关的
   真实信号(情境 B 中滑窗路线吃掉 92% 信号)——flat 模式本质上长期停留在
   该文档明确标注的"受控妥协"区;
2. 无 spread 列,fancy 权重缺失,bin 权重只有样本内项(§3.4);
3. 噪声最大,低幅度事件(单拷贝增益/杂合缺失的边缘)最易漏检。

**适用**:完全没有 normal 可用时的兜底(WGS 且样本自身偏差小者表现相对好);
结果解释需最保守。

## 一句话选型

> 有 matched normal 且要体细胞事件 → `paired`;有一批正常样本或要队列标尺 →
> `pooled`;什么都没有 → `flat`(兜底)。拿不准就留 `auto`,让 node.py 按
> design 决定。

## 配置

```json
"Params": {
    "cnvkit": {
        "enabled": true,
        "method": "wgs",
        "mode": "auto",
        "control_samples": []
    }
}
```

- `method`:CNVkit 的测序方法(`wgs`/`hybrid`/`amplicon`),node.py 会按 meta
  `assay` 列自动推导并覆盖(exome→hybrid、wgs→wgs、panel→amplicon);
- `mode`:`auto`/`paired`/`pooled`/`flat`,决策逻辑见上;
- `control_samples`:pooled 模式的 normal panel 样本名列表(样本名,不是路径;
  样本必须在 meta 中,规则按 `indir/<sid>/<sid><bam_suffix>.bam` 解析其 BAM);
- 配对映射不需要手填:node.py 从 meta design 列推导
  `Params.workflow.paired_design`(`[{"control": normal, "experimental": tumor}]`),
  subworkflow 原样传给模块。

混合方法(cbs 分割、diagram/scatter 图)三种模式一致,由 batch 统一携带。

## 输出

```
{outdir}/
├── reference/reference.cnn          # pooled 模式的共享参考(其余模式无)
└── {sample_id}/
    ├── {sample_id}.cnr              # bin 级 log2 + weight
    ├── {sample_id}.cns              # 片段级 log2
    ├── {sample_id}-scatter.png      # 全基因组散点图
    ├── {sample_id}-diagram.pdf      # 染色体示意图
    └── *.cnn                        # 覆盖度中间表
```

paired 模式下 `{sample_id}` 是配对中的 tumor(normal 不产出 CNV 结果);
pooled/flat 模式下是全部样本。日志按
[modules/modules.md](../modules.md) 日志规范:
单样本规则(paired、batch 内 per-sample 无)落 `logdir/{sample_id}/`,
多样本规则(reference、batch)落 `logdir/cnvkit/`。

## 规则清单

| 规则 | 模式 | 粒度 | 作用 |
|---|---|---|---|
| `cnvkit_reference` | pooled | 组级 | control_samples 池化建 reference.cnn |
| `cnvkit_batch` | pooled/flat | 组级(一条命令全部样本) | 对共享参考/flat 调用 |
| `cnvkit_paired` | paired | 单样本(每对一个 job) | tumor + 自己的 normal 单独调用 |
| `cnvkit_result` | 全部 | — | 面板规则,聚合各模式产物供 node.py 注册 |
