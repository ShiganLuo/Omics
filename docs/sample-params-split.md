# sample_params 拆分(样本参数 vs 流程控制参数)

状态:已实现(2026-11),硬切换、无新旧兼容层。本文记录拆分设计、字段归属、
comparisons 命名归一与迁移影响面。

## 动机

`Params.workflow` 混着两类语义:

1. **流程总控制参数**:aligner / trimmer / snv_caller / counters / input_mode /
   quantification_method / indir / outdir / logdir 等——描述"怎么跑",与样本无关。
2. **样本派生参数**:samples / paired_samples / ip_samples / paired_design /
   comparisons 等——描述"跑什么",**跟随样本数与实验设计变化**,不该住在工具
   参数(Params)里。comparisons 尤其明显:一个 2×2 设计 4 个对比,换个样本集
   整个字段就变。

另外,对比定义历史上有三个同义字段名(`Params.DESeq2.group_pairs` /
`Params.sv.comparisons` / `Params.function.group_pairs`),本次统一为
`sample_params.comparisons` 一个名字、一个形状。

## 新结构

```
{
  "Params": {
    "workflow": { indir, outdir, logdir, meta, dna_indir,
                  aligner, aligner_TEtranscripts, trimmer, snv_caller,
                  counters, input_mode, quantification_method, ... },   # 流程控制(不动)
    "sv": {...}, "DESeq2": {...}, ...                                  # 工具参数(不动)
  },
  "sample_params": {                    # 一级字段:样本派生参数
    "samples": [...], "paired_samples": [...], "single_samples": [...],
    "raw_files": [...],
    "ip_samples": [...], "input_samples": [...],
    "treated_ip_samples": [...], "treated_input_samples": [...],
    "sample_ip_input_map": {...},
    "genome_paired_samples": {...}, "genome_single_samples": {...},
    "dna_paired_samples": [...], "dna_single_samples": [...],
    "paired_design": [...], "pair_strains": {...}, "control_samples": [...],
    "dr_pairs": {...}, "drclass_samples": {...}, "drclass_strain_vcfs": {...},
    "sample_bams": {...}, "populations": {...}, "sample_data": "...",
    "expdesign": "...", "cellranger_dirs": [...], "cellranger_input_dict": {...},
    "comparisons": {
      "<organism>": {
        "<contrast_name>": {
          "control_group_name": "...", "experimental_group_name": "...",
          "control_samples": [...], "experimental_samples": [...]
        }
      }
    }
  },
  ...
}
```

### 字段归属判据

- **sample_params**:值随样本集/实验设计变化(样本名列表、样本角色、配对关系、
  对比、逐样本文件映射),包括 `raw_files`(逐样本收集的输入文件清单)。
- **Params.workflow**:怎么跑(alignment/trim/计数器/输入模式等)与输入输出位置
  (indir/outdir/logdir/meta/dna_indir、spatial 的 input_h5ad/visium_h5/spatial_dir)。
  `genome_pairs`(CoCulture 的双基因组设置)是静态流程设置,不随样本变,
  留在 Params.workflow。
- 模块内部参数(如 `Params.StringTie.sample_groups`、`Params.scanpy.*`)是
  工具适配层的派生值,不算样本注册表,保持原位。

### comparisons 归一(命名 + 形状)

历史上对比定义散在三处、形状各异,均源自同一份 meta group_pairs:

| 旧字段(已废除) | 旧形状 | 消费方 |
|---|---|---|
| `Params.DESeq2.group_pairs` | `{organism: {对比名: {...}}}` | DESeq2 / DR RNA 链 |
| `Params.sv.comparisons` | `{对比名: {...}}`(无 organism) | modules/sv |
| `Params.function.group_pairs` | 各流程不一 | modules/function |

统一为 `sample_params.comparisons`,**organism 嵌套**(多物种无名冲突、无损),
值形状与旧 DESeq2.group_pairs 逐字段一致(control_group_name /
experimental_group_name / control_samples / experimental_samples)。

来源优先级:**config 手填 `sample_params.comparisons` > `run.py --comparison`
文件 > meta design 推导**(手填优先、推导只补缺)。`--comparison` JSON 的
dict 形式按 comparison_id 平铺,进入 config 时按样本 organism 归组。

各模块的对比读取(`comparisons` 模块配置键,统一命名):

| 模块 | 形状 | 接线适配 |
|---|---|---|
| DESeq2 / function(polygenomes) | `{organism: {对比名: {...}}}` | 原样直通 |
| sv(PacVar / Fiberseq) | `{对比名: {...}}` 平铺 | subworkflow 嵌套展平 |
| function 单基因组(Fiberseq) | 平铺 + **组名排序定角色**(匹配 diff_accessibility 约定,旧 Params.function.group_pairs 语义) | Fiberseq.smk `_fn_sorted_comparisons()` |

## 迁移影响面(已完成)

- **写**:node.py 全部 run* 的样本字段写点 → `sample_params`;对比写块统一
  写 `sample_params.comparisons`(runRNAseq / runPacVar / runFiberseq / runDR)。
- **读**:subworkflow/*.smk 全部样本键读点 → `config.get("sample_params", {})`;
  comparisons 接线见上表。模块内的裸键读取(模块配置键名)不变。
- **config 模板**(config/*.json):样本字段移入顶层 `sample_params`;
  旧 comparisons 示例块移除(消费方模板补 `"comparisons": {}` 占位)。
- **schema**(*.schema.json):样本字段声明移入 `sample_params` 节;
  旧 comparisons 声明移除。
- **run.py**:`--comparison` 说明与 `_expand_wildcard_rule_targets` 读点同步。

### 破坏性变更(有意为之,不做兼容)

1. 手填对比的旧位置(`Params.DESeq2.group_pairs` / `Params.sv.comparisons` /
   `Params.function.group_pairs`)失效,一律改填 `sample_params.comparisons`
   (organism 嵌套形状)。
2. 手填样本字段的旧位置(`Params.workflow.<样本键>`)失效,改填
   `sample_params.<样本键>`(如 `sample_params.control_samples`、
   `sample_params.expdesign`)。
3. raw.json 中 `Params.workflow` 不再出现样本键;`sample_params` 是唯一样本
   注册表。

## 测试与验证

- `run.py --test shell -w <全部流程>` 全量回执生成;与重构前回执 diff(路径
  归一化后逐字节一致)即证明传参语义未变(README「回执对照」纪律)。
- `docs/sample-params-split.md` 本文即字段归属与形状的唯一规范来源。
