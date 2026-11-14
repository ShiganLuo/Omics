# openms/raw2mzml — 质谱 RAW → mzML 格式转换子模块

把 Thermo `.raw` 质谱原始文件转换为开放格式 `.mzML`,供 OpenMS 下游
搜谱/定量子模块消费;输入已是 mzML 时直接符号链接。纯工具封装模块。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `samples` (config) | 样本 ID 列表 |
| `raw_files` (config,顶层) | 可选:与 `samples` 一一对应的原始文件完整路径列表(优先于 indir 探测,取 `os.path.realpath`) |
| `indir` (config) | 原始文件目录;`raw_files` 缺席时按序探测 `{indir}/{sample_id}/{sample_id}.raw` / `.RAW` / `.raw.gz` / `.RAW.gz` / `.mzML`,取第一个存在的 |
| `outdir` / `logdir` (config) | 输出/日志目录 |

环境共享父目录:`conda: "../openms.yaml"`、容器 `sif("../openms.yaml", ...)`
指向 `modules/openms/openms.yaml`。

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.mzML     # 转换结果(或对已是 mzML 输入的符号链接)
{outdir}/{sample_id}/raw2mzml.sh          # 回执
```

日志:`{logdir}/{sample_id}/raw2mzml.log`。`raw2mzml_result` 为聚合规则
(expand 全部样本 mzML)。

### 1.3 参数

config 位于 `Params.raw_to_mzml` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.raw_to_mzml.converter` | 引擎可执行 | `thermorawfileparser`(示例 config: `ThermoRawFileParser`) | 转换器 |
| `Params.raw_to_mzml.mode` | — | `thermorawfileparser` | 决定命令形态:`msconvert` / `thermorawfileparser` / 其他(原样 `converter in out`) |
| `Params.raw_to_mzml.args` | 原样前置附加 | "" | 透传转换器参数(空白拆分) |
| `Params.raw_to_mzml.peak_picking` | 见说明 | true | thermorawfileparser 模式:True 默认峰提取,False 追加 `-p`(保留 profile);msconvert 模式:True 追加 `--filter "peakPicking true 1-"` |
| `Params.raw_to_mzml.manifest` | — | "" | 仅出现在示例 config,.smk 未使用 |
| (threads) | — | 2 | `rule_res("raw2mzml")` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
# thermorawfileparser 模式(-f 2 = indexed mzML)
ThermoRawFileParser -i sample.raw -b out/sample/sample.mzML -f 2

# msconvert 模式
msconvert sample.raw --mzML --filter "peakPicking true 1-" \
    -o out/sample --outfile sample.mzML
```

注:thermorawfileparser 固定 `-f 2`(indexed mzML)——pwiz 系读取器
(如 CometAdapter 内部 Comet 包装)对非 indexed 变体(-f 1)会告警并
倾倒文件尾。

## 2. 算法原理与生物学意义

质谱厂商原始格式(Thermo `.raw`)封闭且仪器绑定,mzML 是 PSI 开放
标准(XML,含谱图、保留时间、前体与仪器元数据),几乎所有开源搜谱/
定量工具(OpenMS、Skyline、pyOpenMS)都以它为交换格式。本模块只做格式
转换与峰提取(profile → centroid),不改变谱图数据本身;`peak_picking`
决定是否在转换时做 centroid 化——搜谱一般要 centroid,而部分定量/
去卷积工具偏好保留 profile。已是 `.mzML` 的输入仅做符号链接,保持
DAG 依赖语义一致。

## 3. 参考文献

1. OpenMS — https://github.com/OpenMS/OpenMS(https://openms.de)
2. ThermoRawFileParser(compomics)—— https://github.com/compomics/ThermoRawFileParser
3. ProteoWizard / msconvert — http://proteowizard.sourceforge.net
