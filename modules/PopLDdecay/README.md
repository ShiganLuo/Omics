# PopLDdecay — 群体连锁不平衡(LD)衰减分析模块

封装 PopLDdecay:按亚群计算 SNP 两两间连锁不平衡 r² 随距离的衰减曲线,
输出距离分箱的 r² 统计(`.stat.gz`),用于群体遗传学中有效群体大小、
选择信号与重组历史的评估。单规则模块,引擎为 PopLDdecay 二进制。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `input_vcf`(config) | 群体 VCF(bgzip + tabix 索引)。缺省回退 `outdir/variants/filtered/population_snps.vcf.gz`,同时消费 `.tbi` 索引 |
| `populations`(config) | 亚群字典,键为亚群名;每个亚群对应 `outdir/metadata/populations/{population}.txt`(SubPop 文件:每行 `样本ID<TAB>亚群名`) |
| `Params.analysis.ld_decay` | bool 总开关;false 时规则不定义,模块不进 DAG |

### 1.2 输出

```
{outdir}/analysis/ld_decay/{population}.stat.gz        # 每亚群 LD 衰减统计
{outdir}/analysis/ld_decay/PopLDdecay_{population}.sh  # 命令回执
```

(无 `_result` 聚合规则;下游按 `{population}.stat.gz` 消费。)

### 1.3 参数

config 位于 `Params.analysis` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `ld_decay` | — | false | 模块总开关;false 时不产生规则 |
| `ld_max_distance` | `-MaxDist` | 500000 | 参与统计的 SNP 对最大间距(bp) |

`Procedure.PopLDdecay` 可覆盖可执行路径(缺省 `PopLDdecay`)。

直接调用示例(脱离流程):

```bash
PopLDdecay -InVCF population_snps.vcf.gz \
  -OutStat analysis/ld_decay/CEU \
  -SubPop metadata/populations/CEU.txt \
  -MaxDist 500000
```

产出 `analysis/ld_decay/CEU.stat.gz`。

## 2. 算法原理与生物学意义

PopLDdecay 读取 VCF,在每个 `-SubPop` 定义的亚群内计算物理距离在
`-MaxDist` 范围内的 SNP 对的连锁不平衡 r²(基于基因型计数),按距离
分箱聚合输出 `Dist<TAB>r²` 统计表。生物学上,LD 随距离衰减越快,说明
群体历史有效群体越大、重组打断关联越充分;衰减慢常见于瓶颈、近期扩张
或选择清除区域。读法:横轴 Dist(bp),纵轴平均 r²;曲线半衰点(衰减到
初始值一半的距离)是常用比较量。本模块为纯工具封装。

## 3. 参考文献

1. PopLDdecay — https://github.com/BGI-shenzhen/PopLDdecay
2. Zhang C et al. *PopLDdecay: a fast and effective tool for linkage
   disequilibrium decay analysis based on variant call format files*.
   Bioinformatics, 2019. DOI: 10.1093/bioinformatics/bty875
