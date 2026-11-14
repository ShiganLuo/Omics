# TEtranscripts/disambiguate — 跨物种(移植瘤)TE 表达定量子模块

对 Disambiguate 物种拆分后的 BAM 分别跑 TEtranscripts 家族的 TEcount /
TElocal,定量转座元件(TE)表达并合并各样本计数表。conda 环境共享父目录
`../TEtranscripts.yaml`(容器同源);合并脚本共享父目录
`modules/TEtranscripts/bin/{combineTE.py,geneId2Name.py}`(经
`config["ROOT_DIR"]` 拼接)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config,默认 `output/disambiguate`)+ 样本 ID 布局(形态 A) | `{indir}/{sample_id}/{sample_id}.disambiguatedSpecies_{genome}.bam`:Disambiguate 按物种拆分后的 BAM(`{genome}` 为 genome wildcard) |
| `genome_pairs` | 二元组 `[genomeA, genomeB]`(如 `["hg38","mm10"]`);genome wildcard 只允许这两个值,否则抛 `ValueError` |
| `genome.<G>.gtf` | 每基因组基因注释 GTF |
| `genome.<G>.TE_gtf` | TE 注释 GTF(TEcount 的 `--TE`) |
| `genome.<G>.TEind` | TE 索引(TElocal 的 `--TE`) |
| `samples` | 样本 ID 列表(供合并规则聚合) |

### 1.2 输出

```
{outdir}/TEcount/{genome}/{sample_id}.TEcount.cntTable   # 单样本 TEcount 计数表
{outdir}/TEcount/{genome}/all_TEcount.tsv                 # 全样本合并(ID)
{outdir}/TEcount/{genome}/all_TEcount_name.tsv            # 全样本合并(基因名)
{outdir}/TElocal/{genome}/{sample_id}.TElocal.cntTable    # 单样本 TElocal 计数表
{outdir}/TElocal/{genome}/all_TElocal.tsv / all_TElocal_name.tsv
{logdir}/{sample_id}/{genome}/TEcount.log / TElocal.log
{logdir}/all/TEtranscripts/{genome}_combine_TE*.log
```

收尾规则 `TEtranscripts_result` 聚合 `all_TEcount.tsv` + `all_TElocal.tsv`。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.TEcount` | 可执行文件本身 | `TEcount` | TEcount(TEtranscripts 套件)可执行路径 |
| `Procedure.TElocal` | 可执行文件本身 | `TElocal` | TElocal 可执行路径 |
| `genome.<G>.gtf` | `--GTF` | — | 基因注释 GTF(两工具共用) |
| `genome.<G>.TE_gtf` | `--TE`(TEcount) | — | TE 注释 GTF(家族级定量) |
| `genome.<G>.TEind` | `--TE`(TElocal) | — | TE 位点索引(位点级定量) |
| `genome_pairs` | — | — | 二元组,合法 genome wildcard 白名单 |
| `indir` / `outdir` / `logdir` | — | `output/disambiguate` / `output` / `log` | 输入/输出/日志根目录 |
| `ROOT_DIR` | — | `./` | 项目根目录,定位父目录 bin/ 脚本 |

资源默认:TEcount/TElocal threads=2、mem_mb=32000;combine 规则
threads=1、mem_mb=8000(均经 `rule_res()` 可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
TEcount --sortByPos --format BAM --mode multi \
  -b s1.disambiguatedSpecies_hg38.bam \
  --GTF hg38.gtf --TE hg38_TE.gtf \
  --project s1.TEcount --outdir output/TEcount/hg38

TElocal --sortByPos -b s1.disambiguatedSpecies_hg38.bam \
  --GTF hg38.gtf --TE hg38_TEind \
  --project s1.TElocal        # 产物 s1.TElocal.cntTable 被 mv 到目标路径

python $ROOT_DIR/modules/TEtranscripts/bin/combineTE.py -p TEcount -i output/TEcount/hg38 -o all_TEcount.tsv
python $ROOT_DIR/modules/TEtranscripts/bin/geneId2Name.py -c all_TEcount.tsv -g hg38.gtf -o all_TEcount_name.tsv
```

## 2. 算法原理与生物学意义

TEtranscripts 家族把 RNA-seq reads 归属到转座元件:TEcount 在 TE 家族/
亚家族层面聚合计数(reads 可多重归属,`--mode multi` 多重计入),TElocal
在基因组位点层面计数(用 TEind 索引区分同家族不同插入位点)。二者互补:
家族级稳健但不知位点,位点级精细但计数稀疏。`--sortByPos` 保证输入 BAM
按坐标有序时的确定性计数。生物学意义:移植瘤/跨物种样本先经 Disambiguate
把 reads 按物种拆分,再分物种定量 TE 表达——TE 去抑制是肿瘤基因组不稳定与
应激响应的重要标志;host vs graft 的 TE 表达差异可指示物种特异的调控变化。
结果阅读:cntTable 行为 TE/基因、列为样本;`all_*_name.tsv` 把 geneId 映射
为基因名;下游可做差异 TE 分析(如 DESeq2)。

## 3. 参考文献

1. TEtranscripts / TEcount / TElocal ——
   https://github.com/mhammell-laboratory/TEtranscripts
   (TElocal:https://github.com/mhammell-laboratory/TElocal)
2. Disambiguate —— https://github.com/AstraZeneca-NGS/Disambiguate
   (移植瘤 reads 物种归属;仓库已 archived)
3. Ahdesmäki MJ, Gray SR, Johnson JH, Lai Z. *Disambiguate: an open-source
   application for disambiguating two species in next generation sequencing
   data from grafted samples*. F1000Research 2016, 5:2741.
   DOI: 10.12688/f1000research.10082.1
