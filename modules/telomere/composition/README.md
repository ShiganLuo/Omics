# composition — 端粒组成/变异基序子模块

端粒**长度 + 重复组成**(canonical TTAGGG、变异基序 TVR、ITS)分析,多样本比较 + violin。
引擎:`bin/telomere_profile.py`(自包含,`-s LABEL=BAM` 多样本,`-f` 多格式出图)。

## 输入输出
- 输入(字典模式):`inputs.composition.bam[sample]` = 每样本 fire.sorted.bam。
- 输出:`{outdir}/telomere/composition/`:`telomere_per_read.tsv`、`telomere_summary.tsv`、
  `telomere_composition.tsv`、`telomere_length_violin.{fmt}`(全样本合并图)。

## 参数(`Params.telomere.composition`)
`min_read_len`(默认 5000)、`max_reads`(默认 2e6);多格式 `Params.telomere.formats`。

## 算法
逐读段统计端粒重复:canonical TTAGGG、变异基序(TVR)、ITS 长度;多样本汇总 + 长度分布
violin(全样本一张图)。TVR/组成回答"端粒组成随龄/样本变化"。
