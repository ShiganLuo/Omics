# Fiber-seq 分析流程文档

## 数据特征

- 平台: PacBio Revio
- 化学试剂: SPRQ (BINDINGKIT=103-496-900, SEQUENCINGKIT=103-726-100)
- SMRT Cell: 103-483-800 (Nx)

## 关键发现: Revio SPRQ 不需要 ft predict-m6a

**自 2025 年起，Revio SPRQ 仪器上的 jasmine 自动在测序时生成 m6A 预测。**

验证方法:
```bash
samtools view input.bam | head -1 | grep -o "MM:Z:[^\t]*"
# 输出: MM:Z:A+a.,14,15,2,6,28,...  ← m6A calls 已存在
```

这意味着:
- ❌ 不需要运行 `ft predict-m6a`（会报错 "Model for BINDINGKIT=103-496-900 not available"）
- ✅ 直接从 `ft add-nucleosomes` 开始（使用 BAM 中已有的 m6A calls）
- ✅ 不需要保留 polymerase kinetics tags（IPD/PL）

参考: https://fiberseq.github.io/quick-start.html

## 流程步骤

```
CCS BAM (含 m6A calls from jasmine)
    │
    ▼
ft add-nucleosomes  ← 核小体定位 + MSP 识别
    │
    ▼
ft fire              ← FIRE 调控元件识别
    │
    ▼
ft extract           ← 提取 BED 格式数据 (m6a/nuc/msp/fire)
```

## 工具版本要求

- fibertools-rs >= 0.13.1（编译时需要 `--all-features` 启用 pytorch backend）
- SIF 容器需要: rust, gcc_linux-64, gxx_linux-64, python, pip, torch==2.9.0+cpu, samtools, pbmm2

## SIF 容器构建

```bash
python workflow/Omics/src/common/util/EnvUtil.py apptainer all \
    -y workflow/Omics/modules/fibertools/fibertools.yaml \
    -o /path/to/env/fibertools --force
```

构建关键点:
1. EnvUtil 的 `ENDOFSCRIPT` 必须顶格（无缩进）—— heredoc 关闭标记
2. `micromamba env create` 需要 `< /dev/null` 防止消费 heredoc stdin
3. pip 安装 torch CPU 版（conda 的 libtorch 缺 C++ 头文件）
4. `LIBTORCH_CXX11_ABI=1`（匹配 pip torch 的 ABI）
5. `CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=x86_64-conda-linux-gnu-gcc`
6. `OPENSSL_DIR` + `OPENSSL_NO_VENDOR`（openssl 头文件）
7. LD_LIBRARY_PATH 写入 `/.singularity.d/env/99-fibertools.sh`（Apptainer 自动 source）

## 旧仪器（Sequel II/IIe）处理

旧仪器不支持 jasmine 自动 m6A calling，需要:
1. CCS BAM 必须包含 average kinetics 信息（IPD/PL tags）
2. 运行 `ft predict-m6a`（需要 pytorch backend）
3. 支持的 BINDINGKIT: 101-364-600 等（非 SPRQ）

## 参考文献

- Stergachis et al., 2020, Science (Fiber-seq 原始论文, DOI: 10.1126/science.aaz1646)
- Jha, Bohaczuk et al., 2024, Genome Research (fibertools-rs, DOI: 10.1101/gr.279095.124)
- https://fiberseq.github.io/
- https://github.com/fiberseq/fibertools-rs/issues/71 (BINDINGKIT=103-496-900 问题)
