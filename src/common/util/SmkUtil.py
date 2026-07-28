import os
import yaml
try:
    from ruamel.yaml import YAML
except ImportError:
    YAML = None
import json
import hashlib
import platform

def get_yaml_path(module_name:str,SNAKEFILE_DIR:str)->str:
    """
    获取配置文件路径
    备注：`configfile` 是相对于“执行命令时所在的目录”，而 `include` 是相对于“当前文件所在的目录”
    """
    module_path = os.path.join(SNAKEFILE_DIR ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path


def conda_env_dir(yaml_path: str, frontend: str = "conda"):
    """
    近似 Snakemake 的 conda env hash 计算方式
    """

    # 用 ruamel.yaml 解析，最大程度模拟 Snakemake
    if YAML is not None:
        yaml_ruamel = YAML(typ="rt")
        with open(yaml_path, "r") as f:
            data = yaml_ruamel.load(f)
    else:
        with open(yaml_path) as f:
            data = yaml.safe_load(f)

    # 只提取 channels 和 dependencies 字段
    channels = data.get("channels", [])
    dependencies = data.get("dependencies", [])

    # pip 依赖需要展开
    pip_deps = []
    new_dependencies = []
    for dep in dependencies:
        if isinstance(dep, dict) and "pip" in dep:
            pip_deps.extend(dep["pip"])
        else:
            new_dependencies.append(dep)

    env_dict = {
        "channels": channels,
        "dependencies": new_dependencies,
    }
    if pip_deps:
        env_dict["dependencies"].append({"pip": pip_deps})

    # 只序列化 channels 和 dependencies，缩进2空格，结尾加换行，严格模拟 Snakemake
    if YAML is not None:
        import io
        buf = io.StringIO()
        yaml_ruamel.default_flow_style = False
        yaml_ruamel.indent(mapping=2, sequence=2, offset=0)
        yaml_ruamel.dump(env_dict, buf)
        normalized = buf.getvalue()
        if not normalized.endswith('\n'):
            normalized += '\n'
    else:
        # fallback: json
        normalized = json.dumps(env_dict, sort_keys=True, separators=(",", ":")) + '\n'

    md5 = hashlib.md5(normalized.encode()).hexdigest()
    print(f"Conda env hash for {yaml_path}: {md5}")
    return md5

if __name__ == "__main__":
    conda_env_dir("/data/pub/zhousha/20260411_RNAseq/workflow/RNA-SNP/modules/disambiguate/disambiguate.yaml")
    pass