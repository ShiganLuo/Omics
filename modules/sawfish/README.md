# sawfish — HiFi SV+CNV 联合鉴定模块

PacBio 官方 sawfish(PacificBiosciences/sawfish,bioconda `sawfish=2.3.1`),
对 HiFi 比对 BAM 做 germline SV + CNV 联合鉴定:`discover`(每样本,候选 SV
局部组装 + 分箱深度 + 草稿 CNV 分段)→ `joint-call`(多样本联合基因分型 +
重分段 + SV/CNV 记录合并)→ SVCLAIM 拆分(本模块 bin 引擎)。

## 输出(`{outdir}/`,workflow 侧为 `variation/germline_cnv/`)

```
discover/{sample_id}/            # sawfish discover 输出目录(+ .done 标记)
joint_call/genotyped.sv.vcf.gz   # 合并 VCF(VCF 4.4,SV+CNV 记录)
cnv/sawfish.cnv.vcf.gz(.csi)     # CNV 集:SVCLAIM 含 D 或 SVTYPE=CNV
cnv/sawfish_cnv_summary.tsv      # 记录计数(total/cnv/sv/DEL/DUP/CNV/BND)
cnv/sawfish.sv.vcf.gz(.csi)      # 断点 SV 集(仅 call_sv=true 时产出)
```

## SV 鉴定开关(默认不开启)

sawfish **没有**"只跑 CNV"的命令行开关(`--disable-cnv` 是反向开关;discover
的 SV 候选组装是固定步骤),因此"不开启 SV 鉴定"在两层实现:

1. **计算层** `--fast-cnv-mode`(Params.sawfish.fast_cnv_mode,默认跟随
   `not call_sv`):只保留对 CNV/大变异有用的大断点证据,跳过小片段组装
   (去掉所有插入和 <1kb 的断点缺失)——sawfish 官方的 "CNV/Large variants
   only" 快速模式;
2. **输出层** `Params.sawfish.call_sv`(默认 false):拆分 VCF 只交付
   SVCLAIM 含 D(深度分段支持)的 CNV 记录;断点 SV 记录(SVCLAIM=J/BND)
   默认丢弃,置 true 才写出 `cnv/sawfish.sv.vcf.gz`。

SVCLAIM 语义(VCF 4.4):D=深度分段支持,J=断点组装支持,DJ=两者(拆分时
同时进 CNV 集与 SV 集)。

## 配置(`Params.sawfish`)

| 键 | 默认 | 说明 |
|---|---|---|
| `enabled` | false | workflow 侧门控(node.py 注册);模块本身不读 |
| `call_sv` | false | true 时额外产出断点 SV VCF,并默认关掉 fast_cnv_mode |
| `fast_cnv_mode` | null | null=跟随 `not call_sv`;显式 true/false 优先 |
| `expected_cn` | null | 期望拷贝数 BED(col5=CN;性染色体非整倍体要配),传递 `--expected-cn` |
| `cnv_excluded_regions` | null | CNV 排除区 BED(着丝粒/卫星/常见多态 CNV),传递 `--cnv-excluded-regions` |
| `treat_single_copy_as_haploid` | false | joint-call `--treat-single-copy-as-haploid`(CN=1 区按单倍体分型) |

注意:官方预生成 expected_cn / cnv_excluded_regions 轨道只覆盖人类
hg38/hg19/hs37d5(GRCh38 推荐 `annotation_and_common_cnv.hg38.bed.gz`);
其他物种(如 mulatta)需自制——expected_cn 至少区分性染色体,excluded
可从 assembly gap/着丝粒注释生成(见 sawfish 仓库 `scripts/cnv_excluded_regions`)。

## 资源(实测参考,sawfish 文档)

discover 30x 样本约 1-5 GiB(1-16 线程);joint-call 每样本每线程约 1 GiB,
多样本按 ≥1 GiB/线程 预算。默认 threads 8 / mem_mb 8000(discover)、
16000(joint-call),`rule_res` 可覆盖。

## 直接调用

```
python modules/sawfish/bin/sawfish_split_vcf.py \
    --vcf out/joint_call/genotyped.sv.vcf.gz --out-dir out/cnv [--keep-sv]
```

拆分规则:CNV 集 = SVCLAIM 含 D 或 SVTYPE=CNV;SV 集 = SVCLAIM 含 J 或
SVTYPE=BND 或无 SVCLAIM。bgzip + bcftools index 由引擎内部完成。

## 环境

`env.sawfish` → `/data/16T/Database/env/sawfish/sawfish.sif`(由
`sawfish.def` + `sawfish.yaml` 构建);容器内含 sawfish/bgzip/bcftools/python3。
