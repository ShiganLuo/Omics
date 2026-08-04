# star_3pass / star_3pass_gene 设计说明

`star_3pass` 和 `star_3pass_gene` 的三轮比对均封装在一个 Snakemake 执行规则中。规则内部通过生成带日志和时间戳的 shell script，按顺序执行 pass 1、pass 2、pass 3；subworkflow 只负责 module 配置、`use rule` 和输出编排。

## `star_3pass`

`star_3p_align` 的单个 `run:` 规则执行：

1. 全基因组 end-to-end alignment。
2. 提取 small-RNA reads 后进行 small-RNA reference end-to-end alignment。
3. canonical 与 non-canonical reads 分别进行 genome local alignment。
4. 合并最终 BAM 并建立 index。

## `star_3pass_gene`

`star_3pg_gene_specific` 从已有 `star_3pass` final BAM 开始，在一个 `run:` 规则中执行：

1. 一次性完成全 BAM 到 gene 的 read assignment，并生成 gene-specific FASTQ 和 genomic reference。
2. 对每个 gene reference 建 STAR index。
3. 对每个 gene 顺序执行三轮 local alignment：后续轮次使用前一轮的 unmapped reads（pass1 -> pass2 -> pass3），每轮有独立的 STAR 参数（`passes.pass1` / `passes.pass2` / `passes.pass3`）。
4. 只合并每个 gene 的第三轮 BAM，然后执行 Tailer。
5. 拼接所有 gene 的 tail CSV。

严格模式默认排除 ambiguous read。批量合并 reference 方案保留为未来备选，但目前配置中明确关闭，不参与运行。

注意：dry-run 或规则注册只能证明 DAG 和命令结构，不能证明实际 STAR 结果或文献等价性。