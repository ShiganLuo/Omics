# regions — 端粒区域子模块

从**参考基因组扫描**的端粒坐标(`scan_assembly_telomere` TSV)构建 `telomere_core` /
`subtelomere` 区域 BED(供深度 QC / 区域分析用)。
引擎:`bin/region_build.py telomere --scan-tsv … -o …bed`(自包含)。

## 输入输出
- 输入:`inputs.regions.scan_tsv[sample]`(或 `Params.telomere.regions.scan_tsv`)= 扫描 TSV。
- 输出:`{outdir}/telomere/regions/{sample}/{sample}.telomere.bed`(telomere_core/subtelomere)。

## 参数(`Params.telomere.regions`)
`subtelomere`(默认 100000)、`min_len`(默认 1000)、`scan_tsv`(扫描 TSV 路径)。

## 算法
扫描组装 contig 末端得端粒重复束坐标,转成 telomere_core(重复束)+ subtelomere(内侧窗)
BED;供 length 深度 QC 用参考端粒坐标统计测序深度(见 length §深度 QC)。
