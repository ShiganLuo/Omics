# karrseq — 邻近连接 contact-pair 构建模块(自研脚本组)

从 STAR 比对(`Aligned.out.bam` + `Chimeric.out.sam`)构建 KARRseq 的
contact pairs:提取读对 → 去 PCR 重复 → 4DN pairs 风格导出(pairix
索引),并产出 ligation 断点归一化的 pairs 与距离过滤 BED。三个 bin/
脚本为自研引擎(`get_STAR_reads.py`、`remove_duplicates.py`、
`pairs_to_bed.py`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | STAR 输出目录(上游 config 的 `outdir`):`{indir}/{sample_id}/{sample_id}_Aligned.sortedByCoord.out.bam` 与 `{indir}/{sample_id}/{sample_id}_Chimeric.out.sam` |
| `samples` (config) | 样本 ID 列表 |
| `ROOT_DIR` (config) | 项目根目录,定位 `modules/karrseq/bin/` 脚本 |
| `outdir` / `logdir` (config) | 输出/日志目录 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.txt.gz                    # 排序后的全量 contact pairs(gzip)
{outdir}/{sample_id}/{sample_id}.dedup.txt.gz              # 去重 pairs
{outdir}/{sample_id}/{sample_id}.dedup.bed                 # 距离过滤后的 BED(pairs_to_bed)
{outdir}/{sample_id}/{sample_id}.dedup.pairs.gz            # 4DN pairs(前 7 列,bgzip + pairix 索引)
{outdir}/{sample_id}/{sample_id}.dedup.ligation.pairs.gz   # ligation 断点归一化 pairs(bgzip + pairix 索引)
```

回执:`karrseq_chimeric_to_pairs.sh`、`karrseq_dedup.sh`、
`karrseq_ligation.sh`(各自样本输出目录);日志
`{logdir}/{sample_id}/karrseq_*.log`。pairix 生成的 `.px2` 索引为
附属产物。`karrseq_result` 为聚合规则。

### 1.3 参数

config 位于 `Params.karrseq` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.karrseq.mapq` | 位置参数(get_STAR_reads / pairs_to_bed) | 1 | 读段最低 MAPQ |
| `Params.karrseq.span` | 位置参数(get_STAR_reads) | 0 | span 阈值(嵌合读段跨度过滤) |
| (threads) | `sort --parallel` | 4 | `rule_res` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

自研脚本 CLI(sys.argv 位置参数):

```bash
# 1) 提取 contact 记录(Aligned 模式:读 BAM 流;Chimeric 模式:读 SAM)
samtools view sample_Aligned.sortedByCoord.out.bam \
    | python modules/karrseq/bin/get_STAR_reads.py Aligned 1 0 > tmp1
cat sample_Chimeric.out.sam \
    | python modules/karrseq/bin/get_STAR_reads.py Chimeric 1 0 > tmp2

# 2) 合并排序(键:染色体对/坐标/链)→ gzip
cat tmp1 tmp2 | sort -k2,2 -k4,4 -k3,3n -k5,5n -k10,10n -k11,11n --parallel=4 \
    | gzip -c > sample.txt.gz

# 3) 去重 → BED / pairs
zcat sample.txt.gz | python modules/karrseq/bin/remove_duplicates.py dedup | gzip -c > sample.dedup.txt.gz
zcat sample.dedup.txt.gz | python modules/karrseq/bin/pairs_to_bed.py 1 > sample.dedup.bed
zcat sample.dedup.txt.gz | cut -f1-7 | bgzip -c > sample.dedup.pairs.gz
pairix sample.dedup.pairs.gz

# 4) ligation 断点归一化(awk 平移 + 同染色体方向归一)→ pairix
zcat sample.dedup.txt.gz \
    | awk '{print $1"\t"$2"\t"$3+$10"\t"$4"\t"$5"\t"$6"\t"$7}' \
    | awk '{if($2==$4){if($3>$5){...交换...}else{print $0}}else{print $0}}' \
    | sort -k2,2 -k4,4 -k3,3n -k5,5n --parallel=4 | bgzip -c > sample.dedup.ligation.pairs.gz
pairix sample.dedup.ligation.pairs.gz
```

注:`pairs_to_bed.py` 的距离窗为硬编码 83–100000 bp;`remove_duplicates.py`
第一个位置参数为模式(规则固定传 `dedup`)。

## 2. 算法原理与生物学意义

从数据形态看(STAR 嵌合读段、ligation 断点、pairix 索引的 pairs 文件),
本模块处理邻近连接(3D 基因组)类实验的 contact-pair 构建:

- **get_STAR_reads.py**:STAR 的 Chimeric 输出承载跨连接点嵌合读段
  (两段比对 = 一对互作位点),Aligned 输出承载常规比对;按 MAPQ 与
  span 阈值过滤后转成 contact 记录,合并按「染色体对 → 坐标 → 链」
  排序,保证 pairs 文件确定性;
- **remove_duplicates.py**:PCR/光学重复折叠(同起止坐标读对留一),
  重复接触不重复计数;
- **pairs_to_bed.py**:再按 MAPQ 过滤,只保留基因组距离 83–100 kb 的
  接触(硬编码窗),输出 BED;
- **ligation 规则**:把 mate 起点平移 span(`$3+$10`)后,对同染色体
  对按坐标大小归一方向(小坐标在前),得到与 ligation 断点语义一致的
  pairs;`cut -f1-7` 转 4DN pairs 格式并 pairix 索引,供 Juicebox /
  pairix 生态消费。

生物学意义:contact pairs 量化染色质空间邻近频率,是 3D 基因组
区室/拓扑结构域分析的输入;BED 版本用于快速区间交叠(如与基因/
调控元件注释求交)。模块内未附 KARRseq 方法文献,实验细节请参考
课题组内部文档。

## 3. 参考文献

1. get_STAR_reads.py / remove_duplicates.py / pairs_to_bed.py 为本仓库
   自研脚本;模块文件内未附 KARRseq 公开文献。
2. pairix(4DN pairs 索引工具)—— https://github.com/4dn-dcic/pairix
