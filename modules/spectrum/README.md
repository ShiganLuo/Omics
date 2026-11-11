# spectrum 说明

体细胞突变谱（SBS96）分析模块：VCF → 三核苷酸上下文计数 → 可视化 → 已知签名 refit → de novo NMF 签名提取。

原理、结果影响因素与设计依据见同目录《突变签名分析_原理与影响因素.md》。

## 规则链

```
somatic_spectrum ──> somatic_signature_refit     （NNLS 拟合 COSMIC 已知签名）
        │
        └──────> somatic_signature_nmf           （de novo NMF 提取 + COSMIC 注释）
```

三个规则各自独立；`spectrum_result` 聚合全部产出。

## 输入文件

| 文件 | 来源 | 说明 |
|---|---|---|
| 体细胞 VCF | `sample_somatic_vcf_dict`（node.py 从 tumor/normal 配对自动派生） | 每样本一个，仅读前 6 列 |
| `genome.fasta` | `genome.references.<asm>.fasta` | 三核苷酸上下文提取 |
| `genome.signature_matrix` | `genome.references.<asm>.signature_matrix` | COSMIC 签名矩阵 TSV（SBS 96 通道） |

签名矩阵随 `SigProfilerAssignment` wheel 分发，已解包至
`/data/16T/Database/cosmic_signatures/{GRCh37,GRCh38,mm9,mm10,mm39}/`，
v3.6（与官网当前版本一致），按基因组 build 选对应文件：

- 人类 GRCh38：`/data/16T/Database/cosmic_signatures/GRCh38/COSMIC_v3.6_SBS_GRCh38.txt`
- 小鼠 mm39（GRCm39）：`/data/16T/Database/cosmic_signatures/mm39/COSMIC_v3.6_SBS_mm39.txt`

WES 数据用同目录 `_exome.txt` 变体（外显子机会数归一化版本）。

## 配置（config/Mutation.json）

```json
"Params": {
  "somatic_spectrum": {
    "signature": {
      "enabled": true,
      "mode": "both",
      "min_mutations": 200,
      "min_exposure_frac": 0.02,
      "family_map": null,
      "max_rank": 8,
      "restarts": 20,
      "seed": 42,
      "rank": null
    }
  }
},
"genome": {
  "references": {
    "GRCm39": { "signature_matrix": "/data/16T/Database/cosmic_signatures/mm39/COSMIC_v3.6_SBS_mm39.txt" }
  }
}
```

| 键 | 默认 | 说明 |
|---|---|---|
| `enabled` | false | 签名分析总开关；false 时 node.py 不注册 refit/NMF 产出，规则不进 DAG |
| `mode` | `both` | `refit` / `nmf` / `both`，控制注册哪组产出 |
| `min_mutations` | 200 | 低于此 SNV 数的样本标记 unreliable（采样噪声不可忽略） |
| `min_exposure_frac` | 0.02 | refit 剪枝阈值，低于此曝光比例的签名逐轮剔除后重拟合 |
| `family_map` | null | 自定义家族映射 TSV（`family<TAB>SBS1,SBS2`）；null 用内置表（Clock-like/APOBEC/MMR/POLE_POLD1/HRD/UV 等） |
| `max_rank` | 8 | NMF 扫描的最大秩 |
| `restarts` | 20 | 每个秩的随机重启次数 |
| `seed` | 42 | 随机种子（可复现） |
| `rank` | null | 强制指定秩；null 自动选择（KL 拐点 + cophenetic 下限） |

## 产出

`somatic_spectrum`（组级视图 + 数据基础）：

| 文件 | 内容 |
|---|---|
| `somatic_spectrum_per_sample.csv` | 96 通道 × 样本计数矩阵（refit/NMF 的输入契约） |
| `somatic_spectrum_sample_stats.csv` | 每样本 SNV 计数与丢弃账目（非 SNV / 未知 contig / 无上下文） |
| `somatic_spectrum.csv` | 96 通道 × 分组计数 |
| `somatic_spectrum_stacked_bar.png` / `_heatmap.png` | 组级谱图 |

`somatic_signature_refit`（P1，已知签名比较）：

| 文件 | 内容 |
|---|---|
| `signature_refit_exposures.csv` | 样本 × 签名曝光（突变数） |
| `signature_refit_exposure_fractions.csv` | 行归一化曝光比例 |
| `signature_refit_family_exposures.csv` | 签名家族级聚合比例 |
| `signature_refit_summary.csv` | 每样本拟合统计：n、reliable、用到的签名数、重构余弦、top 签名 |
| `signature_refit_exposure_stacked_bar.png` / `_exposure_heatmap.png` | 曝光图 |

`somatic_signature_nmf`（P2，de novo 提取）：

| 文件 | 内容 |
|---|---|
| `signature_nmf_signatures.csv` | 提取的签名（96 通道 × K） |
| `signature_nmf_exposures.csv` | 样本 × K 曝光 |
| `signature_nmf_rank_diagnostics.csv` / `.png` | 每秩 KL / cophenetic / 稳定性 + 选中的秩 |
| `signature_nmf_signatures.png` | 提取签名的 96 通道剖面 |
| `signature_nmf_cosmic_similarity.csv` / `_best_match.csv` | 与 COSMIC 的余弦矩阵 + 最佳匹配（跨物种时按"跨物种参考"解读） |

## CLI 等价命令

```bash
# 1. 谱计数（somatic_spectrum 规则内部）
python modules/spectrum/bin/spectrum.py \
    -v experiment_vcf_map.tsv -g experiment_group_map.tsv \
    -r GRCm39.primary_assembly.genome.fa -o somatic_spectrum

# 2. refit（somatic_signature_refit 规则内部）
python modules/spectrum/bin/signature_refit.py \
    -c somatic_spectrum_per_sample.csv \
    -s /data/16T/Database/cosmic_signatures/mm39/COSMIC_v3.6_SBS_mm39.txt \
    -o signature_refit -m 200 -e 0.02

# 3. de novo NMF（somatic_signature_nmf 规则内部）
python modules/spectrum/bin/signature_nmf.py \
    -c somatic_spectrum_per_sample.csv \
    -s /data/16T/Database/cosmic_signatures/mm39/COSMIC_v3.6_SBS_mm39.txt \
    -o signature_nmf -k 8 -n 20 -d 42
```

## 跨物种提示

人类 COSMIC 签名拟合小鼠数据（或反之）仅作**跨物种参考**，不构成同源赋值；
小鼠数据建议以 de novo NMF 为主、COSMIC 注释为辅（详见原理文档 §5.6）。
