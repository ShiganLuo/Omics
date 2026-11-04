# gatk_somatic 模块

GATK Mutect2 体细胞 SNV/INDEL 检测(tumor vs matched normal)。从 gatk_prepare
的 markdup BAM 出发,按 meta design 配对逐对 calling,产出每对的原始
(raw/unfiltered)somatic VCF。

## 检测模式

Mutect2 本身只有一种配对 calling 模式,但"基线"有三层可选增强,语义与
CNVkit 模块的三模式对照理解最直观:

| 层 | GATK 参数 | 本模块 config 键 | 性质 |
|---|---|---|---|
| 配对 normal | `-I normal.bam -normal <sid>` | meta design 列(自动) | 逐对参与 calling,体细胞判定的第一证据 |
| 人群 AF 库 | `--germline-resource` | `genome.references.<组装>.germline_resource` | 公共基线,绑定组装版本(gnomAD 等) |
| 自建 normal panel | `--pon` | `Params.gatk.Mutect2.panel_of_normals` | 队列/实验室属性,与本队列 normal 池绑定 |

三层叠加是 GATK 最佳实践:配对 normal 提供该患者胚系的一手证据;
germline-resource 借全人群频率压掉常见多态;PON 压掉本实验系统特有的
批次伪影(捕获偏差、测序错误热点)。三者都可选,缺席时对应 flag 不出现在
命令行。

**没有 flat 兜底**:无 design 配对时本模块不产出任何 VCF(tumor-only
calling 未接线;需要的话用 `Params.gatk.Mutect2` 扩展,但 tumor-only 的
过滤强烈依赖 germline-resource + PON,不建议裸跑)。

## germline_resource 与 panel_of_normals 的区别

两者都用于压制胚系/伪影,但来源与粒度不同:

- `germline_resource`:人群等位基因频率 VCF(须含 AF 字段,如 gnomAD
  af-only_vcf)。作用于 FilterMutectCalls 阶段——人群 AF 高的位点被倾向
  判为胚系(germline filter)。换组装版本必须换对应人群库,所以它是
  参考基因组属性,放 `genome.references.<组装>` 下随 `genome.default`
  切换。
- `panel_of_normals`:自己一批 normal 样本合并出的"出现过位点" VCF
  (上游用 `GenomicsDBImport` + `CreateSomaticPanelOfNormals` 构建,本仓库
  尚无该规则)。PON 中的位点近乎硬过滤,专治 germline-resource 覆盖不到
  的系统性噪声。它描述的是**你的队列/平台**,与组装无关,放 `Params.gatk.Mutect2`。

没有人群库的物种(如小鼠):germline_resource 留 null,配对 normal +
PON 顶上;两者都没有时,过滤完全依赖配对 normal 的证据。

## 配置

```json
"genome": {
    "references": {
        "GRCh38": {
            "germline_resource": "/data/16T/Database/Reference/human/gnomAD/GRCh38/af-only-gnomad.hg38.vcf.gz"
        }
    }
},
"Params": {
    "gatk": {
        "Mutect2": {
            "panel_of_normals": null
        }
    }
}
```

- `germline_resource`:人群库 VCF(bgzip+tabix 索引,含 AF INFO 字段);
  null 时不传。已就位两份:
  - GRCh38:`/data/16T/Database/Reference/human/gnomAD/GRCh38/af-only-gnomad.hg38.vcf.gz`
    (GATK best-practices 桶的 gnomAD af-only 子集,含 `.tbi`,config 已回填)
  - GRCm39:`/data/16T/Database/Reference/mouse/Ensembl/GRCm39/germline_resource.GRCm39.vcf.gz`
    (Ensembl mus_musculus 加工品:contig 改 chr、裁剪 primary、注入 AF=0.001
    占位——Ensembl 小鼠库无群体 AF,此资源是"已知种系位点清单"语义,
    频率分层不可用;同目录 `build_germline_resource.sh` 可重建);
- `panel_of_normals`:PON VCF;null 时不传。构建入口本仓库暂未提供,
  需外部准备;
- BQSR:`genome.known_sites` + `genome.interval` 同时配置时,输入自动
  切到 `.bqsr.bam`;
- hybrid/panel assay:`genome.targets` 非空时自动加 `-L`(与 cnvkit 共用
  同一 targets BED,由 node.py 按 meta assay 推导)。

## 输出

```
{outdir}/
└── {normal}_vs_{tumor}/
    ├── {normal}_vs_{tumor}.vcf.gz    # Mutect2 raw VCF(未 FilterMutectCalls)
    └── somaticMutect2.sh             # 命令回执(无时间戳版)
```

日志:组级规则,落 `logdir/{normal}_{tumor}_mutect2.log`
(见 [modules/modules.md](../../modules.md) 日志规范)。

注意:输出是 **raw VCF**。`FilterMutectCalls` 未作为独立规则接线——VCF
中的 FILTER 列此时全为 `.`;下游做变异筛选(spectrum 频谱、注释)时应按
INFO 字段(如 `AS_FilterStatus`、`TLOD`、`NLOD`)自行阈值过滤,或等
FilterMutectCalls 规则接入后再用 PASS 集。胚系判定信息(population AF /
PON 命中)以 `POPAF`/`PON` 等 INFO 标签写进 raw VCF,germline-resource
与 pon 传入时才有。

## 与 design 的关系

配对来自 meta 的 design 列(ctrl_TAG / exp_TAG),node.py `runMutation()`
解析为 designPairs:每个 pair 一个 job、一对 VCF;同一个 normal 可配多个
tumor(此时该 normal 被排除出胚系 calling 路径)。somatic VCF 路径同时
喂给 spectrum 模块做突变频谱。
