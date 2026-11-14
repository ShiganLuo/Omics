# TEtranscripts — TE 转座子表达定量模块(TEcount + TElocal)

用 TEtranscripts 套件对 RNA-seq BAM 做转座子(TE)表达定量:TEcount 在
TE 亚家族水平计数,TElocal 在 TE 位点水平计数;`bin/combineTE.py` 把各
样本计数表合并成跨样本矩阵。两条支路独立成规则,`TEtranscripts_result`
聚合两份合并表。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 上游 BAM 目录;input 函数按布局拼 `{indir}/{sample_id}/{sample_id}.bam`(形态 A) |
| `samples`(config) | 样本 ID 列表 |
| `genome.gtf`(config) | 基因注释 GTF(TEcount 与 TElocal 共用) |
| `genome.TE_gtf`(config) | TE 注释 GTF(TEcount 用,TE 亚家族注释) |
| `genome.TEind`(config) | TE 位点索引(TElocal 用) |
| `Procedure.TEcount` / `Procedure.TElocal` | 工具可执行路径,缺省取 PATH 中的 `TEcount` / `TElocal` |

### 1.2 输出

```
{outdir}/TEcount/{sample_id}.TEcount.cntTable   # 每样本 TEcount 计数表
{outdir}/TEcount/all_TEcount.tsv                # 跨样本合并矩阵
{outdir}/TElocal/{sample_id}.TElocal.cntTable   # 每样本 TElocal 计数表
{outdir}/TElocal/all_TElocal.tsv                # 跨样本合并矩阵
```

收尾规则 `TEtranscripts_result` 聚合 `all_TEcount.tsv` + `all_TElocal.tsv`。

### 1.3 参数

无 `Params` 配置段;工具参数在规则内固定。线程/内存经
`resources.<rule>.threads/mem_mb` 配置:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| — | `TEcount --sortByPos --format BAM --mode multi` | 固定 | 按坐标排序 BAM、多映射 reads 多重计数模式 |
| — | `TEcount --project` | `{sample_id}.TEcount` | 输出前缀(输出落 `--outdir {outdir}/TEcount`) |
| — | `TElocal --sortByPos` | 固定 | 同上,位点级 |
| — | `TElocal --project` | `{sample_id}.TElocal` | 输出前缀(产 `{project}.cntTable`,规则内 mv 到声明输出) |
| threads | — | TEcount/TElocal=2,combine=1 | `rule_res` 配置 |
| mem_mb | — | TEcount/TElocal=32000,combine=8000 | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
# TEcount(亚家族级)
TEcount --sortByPos --format BAM --mode multi \
  -b S1.bam --GTF genes.gtf --TE TE.gtf \
  --project S1.TEcount --outdir out/TEcount

# TElocal(位点级)
TElocal --sortByPos -b S1.bam --GTF genes.gtf --TE TE_index \
  --project S1.TElocal
mv S1.TElocal.cntTable out/TElocal/S1.TElocal.cntTable

# 合并各样本计数表
python bin/combineTE.py -p TEcount -i out/TEcount -o out/TEcount/all_TEcount.tsv
```

## 2. 算法原理与生物学意义

- **TEcount**:把 reads 按 `--GTF`(基因)+ `--TE`(TE 亚家族)联合注释
  计数;`--mode multi` 下多重比对 reads 在兼容的 gene/TE 之间多重计数
  (而非随机丢弃/唯一归属),适合 TE 这种高同源序列。输出 cntTable 两列:
  `gene/TE` 特征名 + 计数。
- **TElocal**:同样算法,但 TE 注释换为位点级索引(`--TE` 指向
  `genome.TEind`),回答"哪个基因组位置的 TE 被转录"。
- **合并(`bin/combineTE.py`)**:读入各样本 `*.cntTable`,第二列列名改为
  样本 ID,以 `gene/TE` 为键 inner join 成矩阵,写出
  `all_TEcount.tsv` / `all_TElocal.tsv`。
- **生物学意义**:TE 约占哺乳动物基因组 40% 以上,其转录在发育、衰老、
  疾病与基因组不稳定中扮演重要角色;但 TE 序列高度重复,唯一定位的
  reads 极少——multi 模式保留多重比对证据,得到亚家族/位点级的相对
  活性,后续可做差异表达比较(如 TEtranscripts/DESeq2 下游)。
- **结果怎么读**:合并矩阵行 = gene/TE 特征,列 = 样本,值 = 计数;
  `gene/TE` 列形如 `family:class:...` 的 TE 条目与基因 ID 混排
  (combineTE.py 的 `getTE` 用正则 `^[^:]+:[^:]+:[^:]+:[^:]+$` 挑 TE 行)。

## 3. 参考文献

1. TEtranscripts — https://github.com/mhammell-laboratory/TEtranscripts
   (TEcount 实现)。
2. TElocal — https://github.com/mhammell-laboratory/TElocal(位点级定量)。
3. Jin Y, Tam OH, Paniagua E, Hammell M. *TEtranscripts: a package for
   including transposable elements in differential expression analysis of
   RNA-seq datasets*. Bioinformatics 31(22):3593-3599, 2015.
   DOI: 10.1093/bioinformatics/btv422
