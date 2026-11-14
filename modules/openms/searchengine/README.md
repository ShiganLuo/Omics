# openms/searchengine — 蛋白组数据库搜索子模块(Comet/MSGF+/Sage)

OpenMS 适配器三引擎数据库搜索:mzML + decoy 库 → idXML,内置定量化学
(tmt/lfq/dia)驱动的修饰默认值与同位素标签防呆。C 形态子模块,共享父目录
`../openms.yaml` 环境。同目录 `polygenomes/` 为多基因组版子模块(独立
.smk),本 README 覆盖 `searchengine.smk`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `mzML_dir` + 布局参数 | 形态 A 拼接:`{mzML_dir}/{sample_id}/{sample_id}.mzML` |
| `genome.decoy_fasta` | target+decoy 搜索库;未给/不存在时回退 `{decoy_dir}/genome_decoy.fasta`(decoydatabase 模块产出) |
| `samples` / `mzml_files` | 样本列表 |
| `quantification_method` | 定量化学(`tmt`/`lfq`/`dia`),驱动默认修饰集 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_comet.idXML    # CometAdapter 结果
{outdir}/{sample_id}/{sample_id}_msgf.idXML     # MSGFPlusAdapter 结果
{outdir}/{sample_id}/{sample_id}_sage.idXML     # SageAdapter 结果
```

三引擎规则并存,`Params.search_engine.engine`(默认 `comet`)指明首选引擎,
目标由 subworkflow 按此键注册。日志:`{logdir}/{sample_id}/search_engine_*.log`。

### 1.3 参数

三引擎共用形状(config 位于 `Params.search_engine.<engine>`):

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.search_engine.engine` | — | `comet` | 引擎选择(comet/msgf/sage) |
| `*.threads` | `-threads` | `4` | 搜索线程 |
| `*.precursor_mass_tolerance` | `-precursor_mass_tolerance` | `20` | 母离子质量容差(ppm) |
| `*.fragment_mass_tolerance` | `-fragment_mass_tolerance` | `0.02` | 碎片质量容差(Da) |
| `Params.search_engine.comet.fragment_bin_offset` | `-fragment_bin_offset` | `0` | Comet 碎片 bin 偏移 |
| `Params.search_engine.comet.unmatched_action` | `-PeptideIndexing:unmatched_action` | `remove` | 无匹配肽段处理(PeptideIndexer 默认会退出 13,故默认丢弃) |
| `Params.search_engine.msgf.isotope_error_range` | `-isotope_error_range` | `0,1` | MSGF+ 同位素误差范围 |
| `Params.search_engine.comet.fixed_modifications` | `-fixed_modifications`(可重复) | 见下 | 固定修饰 |
| `Params.search_engine.comet.variable_modifications` | `-variable_modifications`(可重复) | 见下 | 可变修饰 |
| `Procedure.comet` / `msgf` / `sage` | — | `CometAdapter` / `MSGFPlusAdapter` / `SageAdapter` | OpenMS 适配器 |

修饰默认值(按 `quantification_method`,显式配置优先):fixed 一律
`Carbamidomethyl (C)`;variable:tmt = `TMTpro (N-term)`、`TMTpro (K)`、
`Oxidation (M)`,lfq/dia = `Oxidation (M)`。两条强制规则(.smk 内置,经
真实 TMT16 数据验证):① 同位素标签修饰必须且永远排在 variable 列表最前
(CometAdapter 会静默忽略 fixed 标签修饰、且标签修饰排后会被丢弃);
② 配了标签修饰但谱图无报告离子时报错(标签/数据错配不许静默)。
环境:`conda: "../openms.yaml"`。

直接调用示例(等价引擎命令):

```bash
CometAdapter -in s1.mzML -out s1_comet.idXML -database genome_decoy.fasta \
  -threads 4 -precursor_mass_tolerance 20 -fragment_mass_tolerance 0.02 \
  -fragment_bin_offset 0 -PeptideIndexing:unmatched_action remove \
  -fixed_modifications "Carbamidomethyl (C)" \
  -variable_modifications "TMTpro (N-term)" -variable_modifications "TMTpro (K)" \
  -variable_modifications "Oxidation (M)"
```

## 2. 算法原理与生物学意义

数据库搜索把谱图与理论碎片谱匹配打分:Comet(相关性打分,快而稳)、
MSGF+ (生成函数打分,统计严谨)、Sage(Rust 现代引擎,快)。OpenMS
适配器统一输出 idXML,供 Percolator/FalseDiscoveryRate 等下游。修饰化学
必须与实验一致:TMT/iTRAQ 标签质量偏移必须作为修饰建模,且 CometAdapter
对 fixed 标签修饰静默失效(模块用 `_normalize_label_mods` 强制归位并排序);
报告离子校验把"标错 lfq/tmt"的配置在搜索步就大声失败。结果怎么读:
idXML 中每 PSM 有得分与 expect 值;target/decoy 注释由 PeptideIndexer 完成,
下游按 1% FDR 过滤。

## 3. 参考文献

1. OpenMS — https://www.openms.de/ ,https://github.com/OpenMS/OpenMS
   (Röst HL et al., Nature Methods 13:741–748, 2016.
   DOI: 10.1038/nmeth.3959)
2. Comet — https://uwpr.github.io/Comet/(Eng JK et al.,
   J Proteome Res 14(12):5243–5251, 2015. DOI: 10.1021/acs.jproteome.5b00644)
3. MS-GF+ — Kim S, Pevzner PA, *MS-GF+ makes progress in a completely
   database-independent peptide identification*. Nature Communications
   5:5277, 2014. DOI: 10.1038/ncomms6277
