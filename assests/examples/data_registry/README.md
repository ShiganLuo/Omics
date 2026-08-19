# data_registry demo

这个示例演示如何把不同格式的 metadata 和分散在多个目录的 FASTQ 文件统一登记成两张表：

- sample_registry.tsv：一行一个样本
- file_registry.tsv：一行一个文件

运行示例：

```bash
python src/download/build_data_registry.py \
  --metadata assests/examples/data_registry/source_a_metadata.tsv assests/examples/data_registry/source_b_metadata.tsv \
  --mapping assests/examples/data_registry/source_a_mapping.json assests/examples/data_registry/source_b_mapping.json \
  --data-root assests/examples/data_registry/files/source_a assests/examples/data_registry/files/source_b \
  --output-dir /tmp/registry_demo \
  --link-dir /tmp/registry_demo_links
```

mapping JSON 约定：

- source_name: 来源名称
- match_key_field: metadata 里用于和扫描文件匹配的列名
- field_map: 来源列名 -> 标准字段
- constants: 给整张表补固定值

标准字段当前包括：
- sample_id
- data_id
- project_id
- assay_type
- organism
- condition
- group
- replicate

如果 metadata 已经直接给了 `fastq_1` / `fastq_2` / `bam` / `cram`，脚本会优先使用这些路径，而不是目录扫描结果。
