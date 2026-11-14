# fragpipe/msfragger — 蛋白组数据库搜索子模块(MSFragger)

MSFragger 数据库搜索:mzML → pepXML/TSV,再经 OpenMS 转换为 idXML
(酶名清洗、spectra_data 注释、PeptideIndexer target/decoy 重注释)供
下游 PSM FDR/定量。C 形态子模块,共享父目录 `../fragpipe.yaml` 环境与
`modules/fragpipe/bin/make_fragger_params.py`;idXML 步骤借用
`../../openms/openms.yaml` 环境(OpenMS 套件)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 布局参数 | 形态 A 拼接:`{indir}/{sample_id}/{sample_id}.mzML`(raw2mzml 产出) |
| `genome.search_fasta` | 搜索蛋白 FASTA(target+decoy,建议用 decoydatabase 模块产出);`_require` 校验存在 |
| `samples` | 样本列表(result 聚合用) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.pepXML            # MSFragger 搜索结果
{outdir}/{sample_id}/{sample_id}.tsv               # MSFragger 表格结果
{outdir}/{sample_id}/{sample_id}_msfragger.idXML   # OpenMS idXML(下游 FDR/定量用)
```

工程细节(见 .smk 注释):mzML 被硬链到 `{outdir}/{sample_id}/` 工作目录
再搜(MSFragger 结果写在输入旁,避免污染 raw2mzml);pepXML 转换走暂存
副本,并把 MSFragger 的 `stricttrypsin` 酶名改写为标准名(OpenMS EnzymeDB
不识别),同时注入 mzML 的 `spectra_data` 路径(IDMapper 要求)。
日志:`{logdir}/{sample_id}/msfragger.log`、`msfragger_idxml.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `data/raw2mzml` / `output` / `logs` | 目录 |
| `samples` | — | `[]` | 样本列表 |
| `genome.search_fasta` | params 文件 `database_name` | 必需 | 搜索 FASTA |
| `Procedure.msfragger` | `msfragger -Xmx...` | `msfragger` | bioconda 启动器(小写) |
| `Params.msfragger.jar` | `java -Xmx.. -jar <jar>` | `""` | 直连 jar 模式(MSFragger-3.x.jar,论文版本,无需 license),优先于启动器 |
| `Params.msfragger.java` / `java_xmx` | `java -Xmx` | `java` / `32G` | jar 模式 JVM |
| `Params.msfragger.license_key` | `--key` | `""` | 启动器模式 license(非空才传) |
| `Params.msfragger.threads` | params 文件 `num_threads` | `16` | 搜索线程 |
| `Params.msfragger.label` | params 生成 `--label` | `tmt16` | 仅 `tmt16` 或 `label-free`,否则 parse 期报错 |
| `Params.msfragger.unmatched_action` | `PeptideIndexer -unmatched_action` | `remove` | 无匹配肽段处理 |
| `Procedure.id_converter` / `peptide_indexer` | `IDFileConverter` / `PeptideIndexer` | OpenMS 工具 | idXML 转换链 |
| 线程/内存 | — | `threads=16`,`mem_mb=32000` | `rule_res("msfragger_run")` 可覆盖 |

环境:`conda: "../fragpipe.yaml"`(搜索)/ `../../openms/openms.yaml`
(idXML);SIF 由 `env.fragpipe` 映射。

直接调用示例(等价引擎命令):

```bash
python modules/fragpipe/bin/make_fragger_params.py --fasta search.fasta --label tmt16
java -Xmx32G -jar MSFragger-3.x.jar fragger.params s1.mzML
IDFileConverter -in s1.pepXML -out s1_msfragger.idXML -out_type idXML
PeptideIndexer -in s1_msfragger.idXML -out s1_msfragger.idXML \
  -fasta search.fasta -unmatched_action remove
```

## 2. 算法原理与生物学意义

MSFragger 是开放式/半开放式搜索引擎:以碎片索引(fragment index)预筛
候选肽段再精细打分,速度比传统引擎快一个数量级,支持 TMT16 标记定量与
label-free。`make_fragger_params.py` 按 `--label` 生成对应参数文件(tmt16:
固定修饰 TMT 标签、开放搜索容差等;label-free:标准酶切)。pepXML → idXML
链把结果送入 OpenMS 生态(PSMFeatureExtractor/Percolator/FalseDiscoveryRate/
Epifany/IDMapper),`PeptideIndexer` 重注释 target/decoy,保证 FDR 一致。
结果怎么读:`.pepXML/.tsv` 逐谱图 PSM(肽段、修饰、expect value);
`_msfragger.idXML` 供下游 FDR 过滤与蛋白推断;TMT 模式下 tsv 含报告离子
通道强度。

## 3. 参考文献

1. MSFragger — https://github.com/Nesvilab/MSFragger ,
   https://msfragger.nesvilab.org/
2. Kong AT et al. *MSFragger: ultrafast and comprehensive peptide
   identification in mass spectrometry-based proteomics*. Nature Methods
   14(5):513–520, 2017. DOI: 10.1038/nmeth.4256
3. OpenMS — https://www.openms.de/(Röst HL et al., Nature Methods
   13:741–748, 2016. DOI: 10.1038/nmeth.3959)
