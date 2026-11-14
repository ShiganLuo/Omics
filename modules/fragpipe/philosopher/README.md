# fragpipe/philosopher — 蛋白鉴定结果验证/过滤模块(Philosopher)

FragPipe 工具链的 Philosopher 步骤:构建含诱饵序列与污染物的搜库
数据库,对上游 MSFragger 产出的 pepXML 做 PeptideProphet 概率验证 +
FDR 过滤(picked protein FDR),导出 PSM/肽段报告。子模块共享父目录
conda 环境(`../fragpipe.yaml`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `genome.search_fasta`(config) | 蛋白序列 FASTA(搜库目标库;必填,`_require` 校验) |
| `indir`(config) | 上游 MSFragger 结果目录(默认 `data/msfragger`);input 函数拼 `{indir}/{sample_id}/{sample_id}.pepXML` |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/database/fragpipe_db.fasta                  # 加诱饵(rev_)+污染物的搜库库
{outdir}/{sample_id}/interact-{sample_id}.pep.xml    # 验证/过滤后的 pepXML
{outdir}/{sample_id}/psm.tsv                         # PSM 级报告
{outdir}/{sample_id}/peptide.tsv                     # 肽段级报告
```

(Philosopher 自产的带日期 `*-decoys-contam-*.fas` 被重命名为
`fragpipe_db.fasta`,保证重跑 mtime 稳定。)收尾规则
`philosopher_result` 聚合全部 `psm.tsv` + 数据库。

### 1.3 参数

config 位于 `Params.philosopher` 下,与 `philosopher filter` CLI 对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `pep_fdr` | `filter --pep` | 0.01 | 肽段级 FDR 阈值 |
| `prot_fdr` | `filter --prot` | 0.01 | 蛋白级 FDR 阈值 |
| `ion_fdr` | `filter --ion` | 0.01 | 离子级 FDR 阈值 |
| `picked` | `filter --picked` | true | picked protein FDR(target-decoy 门控);false 时省略该 flag |
| `Procedure.philosopher` | — | `philosopher` | Philosopher 可执行路径 |
| — | `database --custom ... --contam --prefix rev_ --decoymode 0` | 固定 | 诱饵前缀 `rev_`、附污染物库 |
| — | `peptideprophet --decoy rev_ --nonparam --expectscore --accmass --decoyprobs` | 固定 | 非参数 PeptideProphet 概率模型 |
| threads / mem_mb | — | 1 / 8000(每规则) | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
cd out/database
philosopher workspace --init
philosopher database --custom search.fasta --contam --prefix rev_ --decoymode 0
mv *-decoys-contam-*.fas fragpipe_db.fasta

cd ../S1
philosopher workspace --init
philosopher database --annotate ../database/fragpipe_db.fasta
philosopher peptideprophet --database ../database/fragpipe_db.fasta \
  --decoy rev_ --nonparam --expectscore --accmass --decoyprobs S1.pepXML
philosopher filter --pepxml interact-S1.pep.xml \
  --pep 0.01 --prot 0.01 --ion 0.01 --picked
philosopher report
```

## 2. 算法原理与生物学意义

- **数据库构建**:`database --custom` 把目标蛋白库加上 `rev_` 前缀的
  反序诱饵与常见污染物,组成 target-decoy 库——FDR 估计的基准。
- **PeptideProphet**:对 MSFragger 打分做非参数概率建模
  (`--nonparam --expectscore`),结合精确质量(`--accmass`)与 decoy
  概率(`--decoyprobs`)给每个 PSM 一个"正确鉴定概率"。
- **filter + report**:按 `--pep/--prot/--ion` 三档 FDR 阈值过滤
  (`--picked` 用 picked protein FDR:每个 target 蛋白与其 decoy 竞争,
  蛋白级 FDR 更保守可靠),`report` 导出 `psm.tsv` / `peptide.tsv` 等表。
- **生物学意义**:搜库引擎打分不等于可信度;概率验证 + FDR 控制把
  "多大比例的鉴定是假阳性"钉在阈值内,是定量(如 label-free)前的
  必经步骤。结果读法:`psm.tsv`/`peptide.tsv` 中 Probability/expect
  值高、通过 FDR 过滤的行进入下游定量;`interact-*.pep.xml` 供
  下游(如 MSstats/quantification)继续消费。

## 3. 参考文献

1. da Veiga Leprevost F, Haynes SE, Avtonomov DM, et al. *Philosopher: a
   versatile toolkit for shotgun proteomics data analysis*. Nature Methods
   17(9):869-870, 2020. DOI: 10.1038/s41592-020-0912-y
2. Philosopher — https://github.com/Nesvilab/philosopher
3. PeptideProphet:Nesvizhskii AI, Keller A, Kolker E, Aebersold R. *A
   statistical model for identifying proteins by tandem mass spectrometry*.
   Analytical Chemistry 75(17):4646-4658, 2003. DOI: 10.1021/ac0341261
4. FragPipe — https://fragpipe.nesvilab.org/
