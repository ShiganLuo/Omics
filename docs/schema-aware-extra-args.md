# Schema 感知的 extra_args 类型矫正与非法值校验

## 背景

`run.py` 支持通过 CLI 传额外参数覆盖 workflow config JSON：

```bash
python run.py -w RNAseq --Params.star.outFilterMultimapNmax 3 --counters S1 S2
```

当前 `parse_args()` 的 extra_args 收集逻辑（run.py 第 148-157 行）：

- `--key=v1`       → `extra_args["key"] = "v1"`        （str）
- `--key v1`       → `extra_args["key"] = "v1"`        （str）
- `--key v1 v2 v3` → `extra_args["key"] = ["v1","v2","v3"]` （list）
- `--key`          → `extra_args["key"] = True`         （bool）

**问题**：值类型完全取决于传了几个值，不参考 schema 中的 `type` 定义。导致：
- `--counters S1` → 字符串 `"S1"`，但 config 期望 list
- `--Params.star.outFilterMultimapNmax abc` → 字符串 `"abc"`，但 schema 定义为 int
- 无法发现非法值（如 bool 字段传了 `"maybe"`）

## 设计目标

1. 根据 schema 中的 `type` 字段自动矫正 extra_args 的值类型
2. 对无法矫正的非法值，报错并给出明确提示
3. schema 中未定义的 key 保持原有 `smart_cast` 行为（不阻断）

## 涉及的 schema type

从现有 schema 文件中提取的类型：

| type     | 示例 schema 字段                          | CLI 示例                      |
|----------|------------------------------------------|-------------------------------|
| `str`    | `genome.fasta`, `Procedure.trim_galore`  | `--genome.fasta /path/to.fa`  |
| `int`    | `Params.star.outFilterMultimapNmax`      | `--Params.star.outFilterMultimapNmax 3` |
| `float`  | `Params.function.lfc_cut`                | `--Params.function.lfc_cut 1.5` |
| `bool`   | `Params.function.enabled`                | `--Params.function.enabled true` |
| `list`   | `counters`, `paired_samples`             | `--counters S1 S2 S3` 或 `--counters S1` |
| `array`  | `Params.DESeq2.group_pairs`（带 items）   | 复杂结构，不通过 CLI 覆盖     |
| `dict`   | `Params.scanpy.tissue_samples`           | 复杂结构，不通过 CLI 覆盖     |
| `object` | `genome`（顶层容器）                      | 不直接覆盖                    |
| `null`   | `ROOT_DIR`, `indir`                      | `--indir null`                |

## 改动清单

### 文件 1: `src/common/util/SchemaValidatorUtil.py`

新增两个方法：

#### 方法 1: `get_field_type(dotted_key)`

```python
def get_field_type(self, dotted_key: str) -> Optional[Dict[str, Any]]:
    """根据点号路径查找 schema 中该字段的完整定义。

    Args:
        dotted_key: 点号分隔的 key，如 "counters" 或 "Params.star.outFilterMultimapNmax"

    Returns:
        该字段在 schema 中的定义 dict，找不到返回 None。
    """
```

遍历策略：
- 顶层 key → `self.schema[key]`
- `genome.xxx` → 先试 `self.schema["genome"]["properties"][xxx]`
                  再 fallback `self.schema["genome"][xxx]`
  （兼容两种 schema 结构：RNAseq 的 properties 嵌套 vs CoCulture 的 genome name 嵌套）
- `Params.xxx.yyy` → 递归进入 Params 子树逐层查找
- 其他嵌套 key → 同理逐层下钻，遇到 `"properties"` 就进 properties

#### 方法 2: `cast_extra_args(extra_args)`

```python
def cast_extra_args(self, extra_args: Dict[str, Any]) -> Tuple[Dict[str, Any], List[str]]:
    """根据 schema 类型矫正 extra_args 的值，同时校验非法值。

    Args:
        extra_args: parse_args() 产出的 extra_args dict

    Returns:
        (矫正后的 dict, 错误消息列表)
        错误列表为空表示全部合法。
    """
```

处理逻辑：

| target_type | 矫正规则 | 非法值判定 |
|-------------|----------|------------|
| `list`/`array` | 单值 → `[value]`；已是 list → 不动；递归矫正每个元素（如有 items.type） | 空值且非 nullable → 报错 |
| `int` | `int(value)`，接受 `"42"` 和 `"-1"` | 非数字字符串 → 报错 |
| `float` | `float(value)`，接受 `"1.5"` 和 `"3"` | 非数字字符串 → 报错 |
| `bool` | 匹配 `true/false/yes/no/1/0`（大小写不敏感） | 不匹配 → 报错，提示有效值 |
| `str` | 保持原样 | 无（CLI 传什么都是合法 str） |
| `dict`/`object` | str → `json.loads` | 解析失败 → 报错 |
| `null` | 保持原样 | 非 nullable 且值为空 → 报错 |
| schema 无定义 | 走 `smart_cast` 原有行为 | 不报错 |

### 文件 2: `run.py`

在 `execute_workflows()` 中，extra_args 注入 workflow_config 之前（约第 567 行），插入：

```python
schema_validator = SchemaValidator()
schema_validator._schema_dir = os.path.join(root_dir, "config")
try:
    schema_validator.load_workflow(wf_name)
    casted, errors = schema_validator.cast_extra_args(args.extra_args)
    if errors:
        for err in errors:
            logger.error(f"[{wf_name}] 参数校验失败: {err}")
        raise ValueError(f"extra_args 校验失败，共 {len(errors)} 个错误")
    args.extra_args = casted
except FileNotFoundError:
    pass  # 无 schema 时 fallback 到原有 smart_cast
```

### 不改动的部分

| 组件               | 原因                                         |
|--------------------|----------------------------------------------|
| `parse_args()`     | 仍按原逻辑收集 extra_args，类型矫正交给下游   |
| `smart_cast()`     | 作为 schema 无定义 key 的 fallback             |
| `dict_set_by_path()` | 注入逻辑不变                                |
| `validate()`       | 已有的 required/nullable 检查不变              |
| `get_path_fields()` | 不受影响                                     |
| `generate_test_paths()` | 不受影响                                  |

## 错误信息格式

```
[RNAseq] 参数校验失败: 字段 'Params.star.outFilterMultimapNmax' 期望 int，收到 'abc'
[RNAseq] 参数校验失败: 字段 'counters' 期望 list，但传入空值
[RNAseq] 参数校验失败: 字段 'Params.function.enabled' 期望 bool，收到 'maybe'，有效值: true/false/yes/no/1/0
```

## 数据流

```
CLI
  → parse_args()                    # extra_args: 可能类型不对
  → SchemaValidator.cast_extra_args()  # 按 schema 矫正 + 校验
  → [校验通过] flat_args update + dot_args dict_set_by_path → workflow_config
  → [校验失败] 报错退出，打印所有错误
  → 写入 JSON → snakemake
```

## 边界情况

1. **schema 未定义的 key**：走 `smart_cast`，不报错。允许用户覆盖 config 中 schema 未约束的字段。
2. **array 的 items 是 object**（如 `DESeq2.group_pairs`）：不做自动转换，这类复杂结构用户应通过 JSON 文件传。
3. **bool 合法值**：`true/false/yes/no/1/0`（大小写不敏感）。
4. **int 接受范围**：只接受可解析为整数的字符串（如 `"42"`），`"42.5"` 不接受。
5. **genome 两种 schema 结构**：`get_field_type` 先尝试 `properties` 子树，再 fallback 直接子树。

## 使用示例

```bash
# list 类型：单值也能正确识别为 list
python run.py -w RNAseq -m meta.tsv -o out/ --conda-prefix .conda --counters S1

# int 类型：自动转为 int
python run.py -w RNAseq -m meta.tsv -o out/ --conda-prefix .conda \
    --Params.star.outFilterMultimapNmax 10

# bool 类型：自动转为 bool
python run.py -w scRNAseq -m meta.tsv -o out/ --conda-prefix .conda \
    --Params.cellranger.count.nosecondary true

# 非法值：报错
python run.py -w RNAseq -m meta.tsv -o out/ --conda-prefix .conda \
    --Params.star.outFilterMultimapNmax abc
# 输出: [RNAseq] 参数校验失败: 字段 'Params.star.outFilterMultimapNmax' 期望 int，收到 'abc'
```
