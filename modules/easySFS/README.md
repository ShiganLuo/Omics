# easySFS — 群体位点频率谱(SFS)投影模块

封装 easySFS:按亚群把 VCF 基因型投影为位点频率谱(site frequency
 spectrum),供 dadi/moments 等群体历史推断工具使用。单规则模块
(条件规则:`Params.analysis.history` 为 true 才定义)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `input_vcf`(config) | 群体 VCF(bgzip + `.tbi` 索引)。缺省回退 `outdir/variants/filtered/population_snps.vcf.gz` |
| `populations`(config) | 亚群定义(非空,缺失时模块加载即报错);配套 popmap 文件 `outdir/metadata/population.tsv`(样本 → 亚群) |
| `Params.analysis.history` | bool 总开关;false 时规则不定义 |

### 1.2 输出

```
{outdir}/analysis/history/easySFS/          # easySFS 输出目录(各亚群 SFS 等)
{outdir}/analysis/history/easySFS.done      # 标记文件(DAG 跟踪的声明输出)
{outdir}/analysis/history/easySFS/easySFS_population.sh  # 命令回执
```

(按规范用 `.done` 标记文件跟踪目录型产物;无 `_result` 聚合规则。)

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.analysis.history` | — | false | 模块总开关 |
| `Params.analysis.history_config` | 追加参数 | "" | 透传给 easySFS 的额外 CLI(按空白拆分,如 `--proj 10,12 -a`) |
| `Procedure.easySFS` | — | `easySFS.py` | easySFS 可执行路径 |

直接调用示例(脱离流程):

```bash
easySFS.py -i population_snps.vcf.gz -p population.tsv \
  -o analysis/history/easySFS --proj 10,12 -a
touch analysis/history/easySFS.done
```

## 2. 算法原理与生物学意义

easySFS 读取 VCF 与 popmap,按亚群统计每个位点的衍生等位基因计数,
并把多位点基因型投影(降采样,`--proj`)到固定染色体数,输出每亚群的
一维/多位点 SFS(dadi `.fs` 格式;`-a` 输出 unfolded SFS,需外群定极)。
生物学意义:SFS 携带群体历史信号——近期扩张使低频位点富集、瓶颈使
中频位点减少;把 SFS 喂给 dadi/moments 可拟合有效群体大小变化、分歧
时间与基因流参数。投影是权衡:降采样丢信息但消除测序深度不均的偏差。
结果读法:输出目录中每亚群一个 SFS 文件,元素为各频率档位的位点计数。
本模块为纯工具封装。

## 3. 参考文献

1. easySFS — https://github.com/isaacovercast/easySFS
2. Gutenkunst RN, Williamson SH, Bustamante CD. *Inferring the joint
   demographic history of multiple populations from multidimensional SNP
   frequency data*. PLoS Genetics 5(10):e1000695, 2009.
   DOI: 10.1371/journal.pgen.1000695(SFS 格式的下游消费方 dadi)
