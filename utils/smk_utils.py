import os

def get_yaml_path(module_name:str,SNAKEFILE_DIR:str)->str:
    """
    获取配置文件路径
    备注：`configfile` 是相对于“执行命令时所在的目录”，而 `include` 是相对于“当前文件所在的目录”
    """
    module_path = os.path.join(SNAKEFILE_DIR ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
