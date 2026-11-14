# common — 全模块共享基础设施(common.smk)

`common/common.smk` 不是工具模块,而是被所有 `.smk` 以
`include: "../common/common.smk"`(按层级 `../../` 等价)引入的共享工具层:
日志、资源解析、容器解析、必填校验、蛋白组标记离子检查。它还统一提供
`ROOT_DIR` / `RECEIPT_DIR` 常量与标准库导入。

## 1. 模块输入输出与参数

### 1.1 输入(config 约定)

| config 键 | 默认 | 说明 |
|---|---|---|
| `ROOT_DIR` | `.` | 项目根目录(run.py 注入);用于 `src/` 导入与跨目录脚本定位 |
| `OUTDIR_ROOT` | `output` | 流程总输出目录(绝对路径);回执备份目录 = `{OUTDIR_ROOT}/log/cmd/`(创建一次,文件名自带规则+样本+时间戳) |
| `resources.default` | {} | 全局默认资源(`threads`/`mem_mb` 等) |
| `resources.<rule_name>` | {} | 按规则名覆盖的资源 |
| `env.<rule_name>` / `env.<yaml_stem>` | — | rule/模块级 SIF 容器路径显式映射 |
| `env.env_dir` | — | SIF 目录根(文件系统回退查找用) |

### 1.2 对外接口(函数清单)

| 接口 | 签名 | 说明 |
|---|---|---|
| `setup_logger` | `setup_logger(name, log_file=...)` | 自 `ROOT_DIR/src/common/util/LogUtil.py` 导入的记录器工厂(ROOT_DIR 不对时报带修复指引的 ImportError) |
| `rule_res` | `rule_res(name) -> dict` | 规则资源解析:`resources.default` 与 `resources[name]` 合并(后者覆盖前者);规则内 `threads: rule_res("x").get("threads", 8)` |
| `sif` | `sif(yaml_filename, rule_name="") -> str` | conda YAML → SIF 容器路径解析(见 2.2 四级优先) |
| `_require` | `_require(value, key, tool="")` | 必填 config 输入校验;缺失时抛 `ValueError` 并指名要填的 JSON 键 |
| `check_reporter_ions` | `check_reporter_ions(mzml_path, max_spectra=400) -> (bool, str)` | mzML 同位素标记(TMT/iTRAQ)reporter 离子健全性检查(见 2.3) |
| `ROOT_DIR` / `OUTDIR_ROOT` / `RECEIPT_DIR` | 常量 | 见 1.1;规则回执脚本统一写 `RECEIPT_DIR` |
| 标准库 | `time, shutil, os, sys, tempfile, subprocess, shlex, re` | 随 include 注入,规则可直接用 |

### 1.3 输出

无文件输出;`include` 时创建 `{OUTDIR_ROOT}/log/cmd/` 目录。

## 2. 算法原理与使用约定

### 2.1 rule_res 资源合并

`merged = resources.default ∪ resources.<rule_name>`(按键覆盖,规则无条目
则全继承 default)。典型配置:

```json
"resources": {
  "default": {"threads": 2, "mem_mb": 8000},
  "fiberseq_track": {"threads": 4, "mem_mb": 65000}
}
```

### 2.2 sif 容器解析四级优先(具体优先于通用,显式配置优先于文件系统)

1. `config["env"][<rule_name>]` —— rule 级显式映射;
2. `config["env"][<yaml_stem>]` —— 模块级显式映射(yaml 去扩展名);
3. `<env_dir>/<rule_name>/<yaml_stem>.sif` —— rule 专属容器目录(文件存在才用);
4. `<env_dir>/<module_dir>/<yaml_stem>.sif` —— 模块级回退
   (`module_dir` = `workflow.basedir` 的 basename)。

1/2 均未命中且无 `env_dir` 时抛带修复指引的 `ValueError`。

### 2.3 check_reporter_ions 判定规则(蛋白组标记检查)

对 mzML 全文件蓄水池抽样 `max_spectra`(默认 400)张 MS2 谱(只看文件头会
被无标记的污染物谱误判),解码 m/z 数组(numpress 等无法解码则跳过),
谱内任一峰落在已知 TMT/iTRAQ reporter m/z(单同位素,内置 20 个质量)
±20 ppm 内计命中;**判为标记数据**当且仅当 ≥3 张抽样谱各有 ≥3 个
reporter 质量峰(未标记数据仅有随机质量巧合,~0.1 峰/谱)。返回
`(是否标记, 明细行)`。

## 3. 参考文献

1. Snakemake 官方文档(include / conda / container 机制)——
   https://snakemake.readthedocs.io/
2. 本文件:`modules/common/common.smk`;`setup_logger` 来源
   `src/common/util/LogUtil.py`
