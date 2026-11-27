# length — 端粒长度子模块(Telogator2 / 组装扫描 / 读段密度 / tidk / 深度 QC)

四种端粒长度测量方法并存:① Telogator2 逐染色体臂端粒长度;② 组装 contig
末端端粒重复扫描(推荐小鼠 30–150 kb 端粒,HiFi 读段 15–25 kb 跨不全);
③ 读段级 k-mer 密度全基因组平均;④ tidk 社区工具组装扫描。模块私有脚本
`bin/scan_assembly_telomere.py`、`bin/read_density_telomere.py`、`bin/plot_combined.py`、
`bin/plot_util.py`、`bin/LogUtil.py`。

另有**诊断脚本** ⑤ `bin/telomere_region_depth.py`:统计端粒区比对测序深度、
量化比对恢复率,**证明"区域比对深度法"测端粒长度不可取**(见 §3.4),本身不估长度。

> **方法选择在哪控制**:模块定义全部 rule、不写条件分支;`Params.telomere.method`
> (列表)决定跑哪几种,由 `node.py` 据此注册输出(见 §1.3)——符合"输出决定流程走向"。

---

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_substring` | 形态 A 拼接:`{indir}/{sample_id}/{sample_id}[.{bam_substring}].bam`(+ `.bai`,Telogator2 需要) |
| `assembly_dir` | hifiasm 组装目录(centromere 模块产出):`{assembly_dir}/{sample_id}/assembly/asm.bp.p_ctg.fa`。assembly_scan/tidk 需要;也可用 dict-mode `inputs.assembly_scan.fasta`/`inputs.tidk.fasta` 每样本指定。**均未配置时直接报错(无回退路径)** |
| `samples` | 样本列表 |
| `comparisons` | 对比设计(平铺 `{对比名: {control_group_name, experimental_group_name, control_samples, experimental_samples}}`,subworkflow 从 `sample_params.comparisons` 嵌套展平,见 `docs/sample-params-split.md`)。合并图据此**对照组样本排最左**并标注组名/角色;缺省不排序 |
| tidk 数据库 | `~/.local/share/tidk/tidk_database.csv`;不存在时先跑 `tidk_init`(`tidk build`)建库 |

### 1.2 输出

```
{outdir}/{sample_id}/telogator2/tlens_by_allele.tsv             # 逐等位端粒长度
{outdir}/{sample_id}/telogator2/all_final_alleles.png           # 逐等位端粒长度图
{outdir}/{sample_id}/telogator2/violin_atl.png                  # 端粒长度(ATL)小提琴图
{outdir}/{sample_id}/assembly_scan/{sample_id}_assembly_telomere.tsv        # contig 末端端粒
{outdir}/{sample_id}/assembly_scan/{sample_id}_assembly_telomere_stats.txt  # 汇总
{outdir}/{sample_id}/assembly_scan/{sample_id}_assembly_telomere_length.{fmt}    # 图2 contig 端长度(多格式)
{outdir}/{sample_id}/read_density/{sample_id}_read_telomere.tsv
{outdir}/{sample_id}/read_density/{sample_id}_read_telomere_stats.txt
{outdir}/{sample_id}/read_density/{sample_id}_read_telomere_content.{fmt}         # 图1 竞争/长度分布(多格式)
{outdir}/{sample_id}/tidk/{sample_id}_tidk_telomeres.tsv        # tidk 结果(原生 tsv 改名)
{outdir}/{sample_id}/tidk/{sample_id}_tidk_telomeres.{fmt}       # 图3 tidk 端粒(多格式)
{outdir}/{sample_id}/region_depth/{sample_id}_telomere_region_depth.tsv        # 每区域比对深度(⑤诊断)
{outdir}/{sample_id}/region_depth/{sample_id}_telomere_region_depth_stats.txt  # 恢复率+可取/不可取判定(⑤)
{outdir}/{sample_id}/region_depth/{sample_id}_telomere_region_depth_diagnostic.{fmt}  # 图4 比对深度(多格式)
{outdir}/.tidk_build.done                                       # tidk 建库标记
{outdir}/read_density/combined.{fmt}                            # 图6 全样本合并分布图
{outdir}/combined_methods.{fmt}                                 # 跨方法总览(4 法全指定时)
```

> **图片格式说明(限制)**:`all_final_alleles.png` / `violin_atl.png` 由
> telogator2 内部直接渲染,**无格式开关**(实测 `telogator2 --help` 只有
> `[ALL_FINAL_ALLELES.PNG]` / `[VIOLIN_ATL.PNG]` 两组坐标轴参数,无任何
> 图片格式选项),因此本模块不提供 `formats` 参数、不做多格式展开,输出
> 固定为 png。其他模块的多格式 `formats` 参数不适用于此处(不伪造支持)。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `logs` | 目录 |
| `bam_substring` | 文件名插入段 | `""` | 见输入 |
| `assembly_dir` | — | `""` | 见输入 |
| `Params.telomere.enabled` | — | `false` | 是否启用端粒步骤 |
| `Params.telomere.method` | — | 全部 4 种 | **列表**:要跑的方法(`telogator2`/`assembly_scan`/`read_density`/`tidk`)。空/缺省=全跑;由 `node.py` 据此注册输出,Snakemake 反推只跑所选。选法与模块无关 |
| `Params.telomere.species` | `-t` | `human` | 非 human 时自动找 telogator2 安装内 `non-human/telogator-ref-<species>.fa.gz`,找不到回退默认 |
| `Params.telomere.assembly_scan_length` | `-L` | `50000` | 组装末端扫描长度(bp) |
| `Params.telomere.n_chrom_arms` | `-n` | `40` | 染色体臂数(读段密度法分母) |
| `Params.telomere.formats` | `-F` | `["png"]` | **列表**:出图格式(png/pdf/svg…,可多选)。各分析图据此出多格式 |
| `Params.telomere.threads` | — | 16(rule_res) | Telogator2 线程 |
| `Procedure.telogator2` / `tidk` | `telogator2 -i` / `tidk search` | — | 可执行路径 |
| `Params.telomere.depth_qc.telomere_bed` | `-d` | 无(**必填才跑 ⑤**) | ⑤诊断的端粒区 BED(参考基因组坐标,`telomere_core`/`subtelomere` 四列)。仅配置此键才注册 depth_qc 交付;未配置不产出(**无回退路径**) |
| `Params.telomere.depth_qc.control` | `-c` | `""` | 对照区域(如 `chr1:120000000-120100000`) |
| 固定 | `tidk search --string TTAGGG` | — | 端粒 motif |

脚本 CLI 补充(每选项均有短旗标 + 长名;`-l/--log` 可选写日志文件,缺省只打 stdout):

- `scan_assembly_telomere.py`:`-f/--fasta`、`-s/--sample_name`、`-o/--output_dir`、
  `-L/--scan_length`、`-w/--window_size`、`-e/--step_size`、`-t/--min_tel_frac`、`-l/--log`;
- `read_density_telomere.py`:`-b/--bam`、`-s/--sample_name`、`-o/--output_dir`、
  `-n/--n_chrom_arms`、`-m/--min_read_length`、`-k/--min_tel_kmers`、`-l/--log`;
- `telomere_region_depth.py`(⑤诊断):`-b/--bam`、`-d/--bed`、`-c/--control`、
  `-s/--sample_name`、`-o/--output_dir`、`-m/--min_recovery`、`-l/--log`;
- `plot_combined.py`(多样本合并图):`-m/--method`、`-s/--sample LABEL=TSV`(可重复)、
  `-c/--col`、`-D/--design`(对比设计 JSON,对照组排最左)、`-o/--outdir`、
  `-b/--out_base`、`-F/--format`、`-l/--log`。

环境:`conda: "telomere.yaml"`(模块目录内,pysam/telogator2/tidk)。

直接调用示例(等价引擎命令 + 模块脚本):

```bash
telogator2 -i s1.bam                                  # 非人:加 -t telogator-ref-<species>.fa.gz
python modules/telomere/bin/scan_assembly_telomere.py \
  -f asm.bp.p_ctg.fa -s s1 -o out/s1/assembly_scan -L 50000
python modules/telomere/bin/read_density_telomere.py \
  -b s1.bam -s s1 -o out/s1/read_density -n 40
tidk build && tidk search --string TTAGGG --output s1_tidk_search \
  --dir out/s1/tidk asm.bp.p_ctg.fa
# ⑤ 比对偏差诊断(不估长度;量化纯端粒区比对恢复率,判定区域深度法可取与否)
python modules/telomere/bin/telomere_region_depth.py \
  -b s1.bam -d telomere_regions.bed \
  -c "chr1:120000000-120100000" -s s1 -o out/s1/region_depth
```

---

## 2. 结果解读

### 2.1 表格 / 文本输出怎么读

- `*_stats.txt`(assembly_scan / read_density):给**平均端粒长度**与分布统计
  (Count/Mean/Median/Std/Min/Max/Q25…);read_density 还给 `Telomeric fraction`
  (端粒碱基占比,覆盖无关)与 `Avg telomere per arm`(= total_tel_bp ÷ n_chrom_arms,
  口径见 §3.3)。
- `tlens_by_allele.tsv`(telogator2):**逐等位**(per-allele)端粒长度,每行一个等位;
  端粒等位间差异大,这是它的价值。
- `*_assembly_telomere.tsv`:每 contig 每端的端粒长(位置/长度/motif/支持窗口数)。
- `*_tidk_telomeres.tsv`:tidk 原生输出,组装上搜到的端粒重复段。
- ⑤ `*_telomere_region_depth.tsv/stats.txt`:每区域比对深度 + 恢复率 + "可取/不可取"判定。

方法间互验:②/③/④ 互为验证;小鼠以 ② 为准,样本 ATL 以 ③ 为准。

### 2.2 图像结构与元件含义

> 各图默认 png,可用 `Params.telomere.formats`(列表)同时出多格式(png/pdf/svg…,
> 规范#9)。telogator2 两图是它自带、固定 png(无格式开关,例外)。

- **图 1 · `{s}_read_telomere_content.{fmt}`(read_density,3 面板)**:
  ①逐读段端粒片段长度直方(横=端粒 bp/log,纵=读段数)——长度代理,看中位/高分位;
  ②末端端粒长度直方(仅锚定读段)——端粒长度代理;n 标注锚定读段数;
  ③全基因组端粒占比柱(= total_tel_bp/总读段碱基,覆盖无关)——**含量**。①②是长度、
  ③是含量,二者分开看(含量≠长度)。
- **图 2 · `{s}_assembly_telomere_length.{fmt}`(assembly_scan,2 面板)**:
  ①contig 端端粒长度直方(≥1kb);②逐 contig 端长度散点(5'/3' 不同色/形)。
  **组装/参考**端粒长(非样本个体);看长度分布集中度、p/q 端差异。
- **图 3 · `{s}_tidk_telomeres.{fmt}`(tidk,2 面板)**:①tidk 端粒重复长度/大小分布;
  ②逐 contig 均值柱。社区工具独立验证组装端粒,与图 2 互验。列名自适应解析。
- **图 4 · `{s}_telomere_region_depth_diagnostic.{fmt}`(⑤诊断,2 面板)**:
  ①各区域类平均比对深度柱(telomere_core/subtelomere/对照);②恢复率条形 + 0.5 阈值线。
  **方法学 QC**:看纯端粒深度骤降、恢复率远低于阈值 → 证明区域深度法测长不可取(非生物学结论)。
- **图 5 · `violin_atl.png` / `all_final_alleles.png`(telogator2,自带)**:
  端粒长度(ATL)小提琴图(中位/分布宽窄);逐等位端粒长度(等位间长度差)。
- **图 6 · `read_density/combined.{fmt}`(read_density 多样本合并,2 面板)**:
  ①逐读段 `tel_bp_estimate` 密度直方(横轴 log,`density` 归一——样本读段数不同,
  计数不可直接比);②各样本小提琴+箱线(纵轴 log),中位数标数值。
  配色按对比设计分组;有 `comparisons` 时**对照组样本排最左**,x 轴标签
  标注「组名 · control/experiment」,图例为唯一图级图例(无标题,轴标签带
  方法/列名)。跨样本比较优先看 `tel_fraction`(§3.3),`tel_bp_estimate`
  受测序量影响。

---

## 3. 算法原理

端粒(TTAGGG 重复)长度是细胞衰老与基因组稳定性的经典指标。难点:**纯端粒
(TTAGGG)n 超出 HiFi 读长、且极难唯一比对**——这直接决定哪些方法可信。模块并行
多种方法,并用诊断脚本 ⑤ 摆出证据说明"比对/区域深度法"为何不可取。

### 3.1 各方法原理、输入数据与适用

| # | 方法 | 输入数据 | 原理 | 输出粒度 | 受比对偏差影响 | 适用 / 推荐 |
|---|---|---|---|---|---|---|
| ① | Telogator2 | 比对 BAM + 索引(`.bai`) | 读段种子-延伸匹配端粒/亚端粒参考 | 逐染色体臂、逐等位 | **是**(走参考比对) | 人/端粒较短;逐臂分布 |
| ② | Approach A 组装扫描 | 组装 FASTA(hifiasm) | contig 两端 `-L` 窗口找连续 CCCTAA/TTAGGG 段 | 每 contig 每端 | 否(在组装上) | **小鼠首选**(30–150 kb);给的是**参考/组装**端粒长,非样本 |
| ③ | Approach B 读段密度 | **原始 / unmapped BAM 或 reads(须含全部读段)** | **全读段** telomeric k-mer 总碱基数 ÷ 染色体臂数 | 全基因组平均 | **否**(不比对,数 k-mer) | **无偏 ATL 首选**;数原始读段,含未比对 |
| ④ | Approach C tidk | 组装 FASTA | 社区工具 `tidk search` 在组装上搜端粒重复 | 每 contig | 否(在组装上) | 组装端粒验证 |
| ⑤ | `telomere_region_depth.py` | **比对 BAM** + 端粒区 BED | 统计端粒区**比对**测序深度 + 恢复率 | 每区域 | 是(仅诊断) | **不估长度**,只证明 ①/区域深度法的偏差 |

> **输入数据怎么选(③ 的关键)**:③ 用 pysam 读 BAM,但数的是**读段序列的 k-mer**
> (`fetch(until_eof=True)`,含 `is_unmapped` 但有序列的读段),方法参考无关——
> **但只在输入 BAM 含全部读段时成立**。故 ③ 必须喂**原始/unmapped BAM(如
> `1_raw_data/*.HiFi.bam`)或 reads**;喂**比对+过滤 BAM(如 `3_align_bam/*.sorted.bam`,
> 实测 unmapped=0、纯端粒区深度≈0)会丢掉纯端粒读段,即使方法参考无关也照样低估 ATL**。
> "bam 天然不合适"要精确成:**比对/过滤 BAM 不合适,原始 unmapped BAM 合适**。
>
> **方法 ② 与 ③ 的区别要分清**:② 扫的是**组装 contig 末端**,给的是参考/组装的
> 端粒长度(不反映样本个体端粒);③ 数的是**样本原始读段**,给样本真实 ATL。
> 测"样本端粒长度"应选 ③。

### 3.2 保留多重比对会怎样?

纯端粒读段在多条端粒间高度相似,比对时会**多重比对**(1 条 primary + N 条 secondary,
flag 0x100),每条记录的 `query_sequence` 是**同一条物理读段**。此时:

- **保留 secondary 一起数 k-mer** → 同一条端粒读段被数 N+1 次 → **端粒量虚高(重复计数)**;
- **丢弃 secondary(如 `-F 0x900`)** → 端粒信号丢失 → **低估**;
- 即:比对 BAM 里无论去留多重比对都有问题——**去则漏、留则重**。

而**原始(unmapped)读段不做比对,每条读段只出现一次**,天然没有多重比对/重复计数问题。
这正是 ③ 用原始 reads 作输入能绕开整个问题的原因。**结论:别在比对 BAM 上纠结多重比对
去留;直接用原始/unmapped 读段。**

### 3.3 n_chrom_arms 与两个指标的区别

③ 报两个量——`Telomeric fraction = total_tel_bp / total_read_bp`(端粒碱基占比,
**覆盖无关、稳健**)与 `Avg telomere per arm = total_tel_bp / n_chrom_arms`
(绝对长度,**覆盖依赖**)。

- `n_chrom_arms` 要按物种/倍性设:模块默认 40(小鼠 20 染色体二倍体计数);
  恒河猴二倍体计数 42(20 对常染色体 + XY;Mmul_10 参考的 22 条序列是
  1–20 + X + Y,X/Y 同属一对性染色体,不能按 22×2 计;与 GRCh38=46、
  GRCm39=40 同一换算口径——人参考 24 条序列同样不按 48 计),若按端粒端数
  (每染色体 2 端 ×2 拷贝)则为 84,口径不同绝对值差一倍。
  node.py 会按 `genome.default` 自动填(`GENOME_N_CHROM_ARMS`:小鼠 40、
  恒河猴 42、人 46),显式配置优先;未收录的基因组必须在配置里手动设置。
- `Avg per arm` 还随测序量线性放大(测得越多 total_tel_bp 越大),**不同深度样本间
  比较应看 `Telomeric fraction`**;同深度样本的 `per arm` 才可粗比。
- 故跨样本比端粒长度,优先用 `Telomeric fraction`;要绝对长度需先按覆盖/单拷贝 k-mer 归一。

### 3.4 比对偏差诊断:为什么"区域比对深度法"测端粒长度不可取

直觉上可在端粒区取比对读段、按覆盖深度/重复量推端粒长度。**但纯 (TTAGGG)n 读段
几乎比对不上**——诊断脚本 ⑤ 实测(恒河猴 ovary HiFi,luanchao-21310):

| 区域 | 平均比对深度 | 相对常染色质 |
|---|---|---|
| `telomere_core`(纯 TTAGGG) | **0.00–0.11×** | **~0–0.3%** |
| `subtelomere` | 33.9× | ~80% |
| 常染色质对照 | 42–43× | 100% |

即**纯端粒区比对深度比基因组平均低 ~300 倍,~99.7% 的纯端粒读段未比对/多比对而丢失**。

**为何不可取:**
1. 端粒长度 = 端粒重复的量;长度信息恰在那 99.7% 丢失的读段里 → 用比对深度/重复量
   估计 ATL 会**系统性低估**,且几乎无数据(0.1× 覆盖)。
2. 残留能比对的读段是带**唯一亚端粒锚点的偏倚子集**,不是随机样本,不能外推。
3. 结论:**区域比对深度法不可取**。测 ATL 用 ③(参考无关 k-mer)或 ①(逐臂);
   组装端粒长用 ②/④。脚本 ⑤ 把恢复率打印出来,低于阈值(默认 0.5)即自动判定"不可取"。

**如何判断某数据是否也中招(通用判据)**:用脚本 ⑤ 比较 `telomere_core` 深度 vs
常染色质对照深度,恢复率 = 二者之比;恢复率越低,比对深度法越不可信(本例≈0)。

> **补充(与 motif 分析的差别)**:比对偏差对**基序普查**影响小——普查只需"重复单元"
> 而非均匀深度,带锚点读段虽稀疏但携带长端粒束(数千 TTAGGG/TVR 单元),足以枚举基序谱;
> 但**逐端粒精确定量与长度估计**则完全不可靠。故 motif 结果可留,长度必须走 ③/①。

---

## 4. 参考文献

1. Telogator2 — https://github.com/telomere-toolkit(Telomere toolkit;
   如引用请以官方页面为准)
2. tidk — https://github.com/tolkit/tidk
3. 端粒与衰老综述:López-Otín C et al. *Hallmarks of aging: An expanding
   universe*. Cell 186(2):243–278, 2023. DOI: 10.1016/j.cell.2022.11.001
