# openms/decoydatabase — 蛋白组搜索 decoy 数据库生成子模块(OpenMS DecoyDatabase)

从蛋白 FASTA 生成 target+decoy 搜索数据库,供 PSM 鉴定的 FDR 估计使用。
C 形态子模块,共享父目录 `../openms.yaml` 环境(OpenMS 套件)。
同目录 `polygenomes/` 为多基因组版子模块(独立 .smk),本 README 覆盖
`decoydatabase.smk`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `genome.protein_fasta` | 蛋白 FASTA(`-in`);规则 input 函数 `_require` 校验存在。模块内置核酸 FASTA 防呆:序列字母表仅 ACGTN 时直接报错(核酸库搜蛋白会静默得到 0 PSM) |

### 1.2 输出

```
{outdir}/genome_decoy.fasta     # target+decoy 搜索库
```

日志:`{logdir}/decoy_database/decoy_database.log`;回执脚本
`{outdir}/decoy_database.sh`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `data/fasta` / `output` / `logs` | 目录 |
| `genome.protein_fasta` | `-in` | 必需 | 蛋白 FASTA |
| `Procedure.openms` | — | `DecoyDatabase` | OpenMS DecoyDatabase 可执行路径 |
| `Params.decoy_database.decoy_string` | `-decoy_string` | `DECOY_` | decoy 序列名前缀/后缀标记 |
| `Params.decoy_database.decoy_string_position` | `-decoy_string_position` | `prefix` | 标记位置 |
| `Params.decoy_database.method` | `-method` | `shuffle` | decoy 生成方法(shuffle 等) |
| `Params.decoy_database.shuffle_max_attempts` | `-shuffle_max_attempts` | `30` | 打乱最大尝试次数 |
| `Params.decoy_database.shuffle_sequence_identity_threshold` | `-shuffle_sequence_identity_threshold` | `0.5` | 打乱后序列一致度上限 |

环境:`conda: "../openms.yaml"`。

直接调用示例(等价引擎命令):

```bash
DecoyDatabase -in proteome.fasta -out genome_decoy.fasta \
  -decoy_string DECOY_ -decoy_string_position prefix -method shuffle \
  -shuffle_max_attempts 30 -shuffle_sequence_identity_threshold 0.5
```

## 2. 算法原理与生物学意义

纯工具封装:OpenMS DecoyDatabase 按 `method`(默认 shuffle:打乱每条蛋白的
残基顺序,保留组成)生成 decoy 序列,加 `decoy_string` 标记后与 target
合并。数据库搜索中 decoy 命中数估计假阳性率,是 PSM/蛋白 FDR(如 1%)过滤
的统计基础;`shuffle_sequence_identity_threshold` 防止打乱后仍与天然序列
过同源导致 decoy 失真。结果怎么读:`genome_decoy.fasta` 中带 `DECOY_` 前缀
的条目即 decoy;下游 searchengine/psmfdr 按前缀区分 target/decoy。

## 3. 参考文献

1. OpenMS — https://www.openms.de/ ,https://github.com/OpenMS/OpenMS
2. Röst HL et al. *OpenMS: a flexible open-source software platform for
   mass spectrometry data analysis*. Nature Methods 13(9):741–748, 2016.
   DOI: 10.1038/nmeth.3959
