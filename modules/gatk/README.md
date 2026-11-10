# gatk说明

## gatk_prepare（gatk_prepare.smk）

规则链：gatk_index（参考索引三件套）→ addReadsGroup（RG）→ MarkDuplicates（排序+标重 BAM，含 Markdup-metrics）→ collectHsMetrics（捕获 QC）。

### collectHsMetrics

对 analysis-ready 的 `{sample}.sorted_markdup.bam`（与变异调用完全相同的文件）按捕获区间运行 GATK CollectHsMetrics，报告 on-target 深度、覆盖均一性（PCT_TARGET_BASES_*X）、fold enrichment、off-bait 率等杂交捕获 QC 指标。

- 输入：markdup BAM/BAI、`genome.targets`（捕获 BED 或 Picard interval_list）、fasta/fai/dict（配置了 `fai_index`+`dict_index` 用外部三件套，否则用 gatk_index 规则的产物）
- BED 经 BedToIntervalList 转为 interval_list（bait=target：区间由实际覆盖推导，无独立 bait 设计）；`.interval_list` 输入直接 staged
- 输出：`{outdir}/{sample}/{sample}.hs_metrics.txt` 与 `{sample}.targets.interval_list`
- 缺 `genome.targets`（如 WGS）时不注册输出、规则不进 DAG；被请求而未配置则建图期 ValueError
- 容器：gatk.yaml → gatk4.sif（sif("gatk.yaml", "collectHsMetrics")，无 per-rule 容器时回落模块级）
- 注意：GATK 4.6.2.0 的 CollectHsMetrics 不支持 `--THREADS`，单线程运行，资源键 `collectHsMetrics`（threads 4 / mem_mb 16000，见 config/Mutation.json）

gatk只适用于NGS数据，不适用于三代测序数据

