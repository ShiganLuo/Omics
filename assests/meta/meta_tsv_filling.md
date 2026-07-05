# meta.tsv 填充通用指南

从 run.tsv 和 fastq 实际路径生成 meta.tsv 的标准流程。

## meta.tsv 标准列定义

| 列名 | 说明 | 来源 |
|---|---|---|
| `sample_id` | 样本唯一标识，用于流程输出目录命名 | 从 run.tsv "Run title" 列提取，若无则用 Accession |
| `data_id` | 数据库登录号（CRR/SRR/ERR 等） | 从 run.tsv "Accession" 列提取 |
| `design` | 实验分组（去除非重复后缀） | 从 sample_id 推断：去掉 `-1`/`-2`/`-3`、`_Rep1`/`_Rep2` 等重复后缀 |
| `fastq_1` | R1 绝对路径（单端则为唯一定向文件） | 实际 fastq 文件路径 |
| `fastq_2` | R2 绝对路径（单端留空） | 实际 fastq 文件路径，单端测序留空 |
| `workflow` | 工作流名称 | 根据数据类型和 parametes.txt 判断 |

## 填充流程

### Step 1: 读取 run.tsv 并解析列结构

run.tsv 的列结构可能因数据来源（GSA/NCBI/DDBJ）而异。常见列映射：

```
GSA 格式: ID / Accession / Run title / ... / Read filename 1 / Read filename 2 / ...
NCBI SRA: Run / ReleaseDate / LoadDate / spots / bases / spots_with_mates / avgLength / ...
```

关键提取字段：
- **样本名**: GSA 用 "Run title"，NCBI SRA 用 "Sample Name" 或自己构造
- **登录号**: GSA 用 "Accession"（CRR 开头），NCBI 用 "Run"（SRR 开头）
- **文件名**: GSA 用 "Read filename 1"（含文件大小），NCBI 用 "download_path"
- **是否 paired**: 检查 "Read filename 2" 列是否非空

### Step 2: 确认 fastq 文件实际位置

```bash
# 列出所有 fastq 文件
find /path/to/fastq_dir -name "*.fastq.gz" -o -name "*.fq.gz" | sort

# 统计数量与 run.tsv 行数对比
ls /path/to/fastq_dir/*.fastq.gz | wc -l
```

文件命名约定：
- `{accession}.fastq.gz` — GSA 下载，文件名即 Accession
- `{accession}_1.fastq.gz` / `{accession}_2.fastq.gz` — paired-end
- `{accession}.single.fq.gz` — 项目内部单端命名

### Step 3: 推断 design 分组

design 列的值 = sample_id 去除重复编号/生物学重复后缀。规则：

1. **末尾数字后缀**: `-1` / `-2` / `-3` → 去掉（`GV-1` → `GV`）
2. **Rep 后缀**: `_Rep1` / `_Rep2` → 去掉（`mA_A_0_4_Rep1` → `mA_A_0_4`）
3. **保留有意义的数字**: `mESC-1ng` 中的 `1ng` 是浓度，不去掉
4. **保守原则**: 不确定时保留原样，让用户确认

常见 pattern：
```
# 发育阶段样本: Stage-Rep → Stage
GV-1, GV-2, GV-3 → GV
MII-1, MII-2 → MII
4C-1, 4C-2, 4C-3 → 4C

# 细胞系+处理量: CellType-Amount-Rep → CellType-Amount
mESC-10ng-1, mESC-10ng-2 → mESC-10ng
MEF-100cell-1 → MEF-100cell

# 下划线分隔 + Rep: Remove _RepN suffix
mA_A_0_4_Rep1, mA_A_0_4_Rep2 → mA_A_0_4
```

### Step 4: 确定 workflow

根据数据类型和 parametes.txt / 实验描述判断：

| 数据特征 | workflow |
|---|---|
| 单端 + UMI 提取 + tRNA 分析 | `tRNAseq` |
| 双端 RNA-seq（mRNA） | `RNAseq` |
| 双端 small RNA / ncRNA | `ncRNAseq` |
| 双端 ChIP-seq / CUT&Tag | `PeakCalling` |
| 双端 WGS / Exome | `Mutation` |
| PacBio HiFi 长读长 | `PacVar` |
| CLIP / iCLIP / eCLIP | `CLIP` |
| MeRIP / m6A-seq | `MERIP` |

### Step 5: 组装并验证

1. 用脚本或手动组装 TSV（tab 分隔）
2. **验证所有 fastq 路径存在**：
   ```bash
   awk -F'\t' 'NR>1 && $4!="" {if(system("test -f "$4)!=0) print "MISSING: "$4}' meta.tsv
   ```
3. **验证 design 分组合理性**：
   ```bash
   awk -F'\t' 'NR>1{d[$3]++} END{for(k in d) print k"\t"d[k]}' meta.tsv | sort
   ```
4. **验证无重复 sample_id**：
   ```bash
   awk -F'\t' 'NR>1{print $1}' meta.tsv | sort | uniq -d
   ```

## 常见 Pitfalls

### 1. fastq_2 列：单端留空，不是填 "NA" 或 "none"
Snakemake 配置解析时，空字符串代表单端。填 "NA" 会导致路径拼接错误。

### 2. design 分组过度聚合
`mESC-1ng-1` 和 `mESC-1ngRNA-1` 的 design 不同（`mESC-1ng` vs `mESC-1ngRNA`），不要错误合并。

### 3. sample_id 必须唯一且不含特殊字符
避免空格、斜杠、冒号。下划线和连字符安全。

### 4. fastq 路径必须是绝对路径
Snakemake 的工作目录可能与 meta.tsv 所在目录不同，相对路径会导致找不到文件。

### 5. Run title 可能为空
GSA 的某些记录没有 Run title，此时用 Accession 作为 sample_id。

### 6. paired-end 的 fastq_1/fastq_2 对应关系
同一 sample_id 的 `_1.fq.gz` 和 `_2.fq.gz` 必须成对出现，路径在 fastq_1 和 fastq_2 两列。

## 输出模板

```tsv
sample_id	data_id	design	fastq_1	fastq_2	workflow
Sample1	SRR000001	GroupA	/path/to/SRR000001_1.fq.gz	/path/to/SRR000001_2.fq.gz	RNAseq
Sample2	SRR000002	GroupA	/path/to/SRR000002_1.fq.gz	/path/to/SRR000002_2.fq.gz	RNAseq
```
