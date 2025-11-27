from itertools import combinations
from statannotations.Annotator import Annotator
import seaborn as sns
import matplotlib.pyplot as plt

# 假设 df 是您的 DataFrame

# --- 步骤 1: 确定所有唯一的组别 ---
unique_cell_types = df['Cell_Type'].unique().tolist()
# 示例: ['senescent cells', 'tumor cells', 'another cell type']

# --- 步骤 2: 自动生成所有两两比对的组合 ---
box_pairs = list(combinations(unique_cell_types, 2))
# combinations(unique_cell_types, 2) 会生成所有不重复的 (组A, 组B) 组合

# --- 步骤 3: 绘制图表并添加差异线 ---
plt.figure(figsize=(10, 8))
ax = sns.violinplot(
    x='Cell_Type',
    y='LV1',
    data=df,
    palette='Set2',
    inner='box',
    linewidth=1.5
)

# 初始化 Annotator，使用自动生成的 box_pairs
annotator = Annotator(
    ax,
    box_pairs,
    data=df,
    x='Cell_Type',
    y='LV1',
    order=unique_cell_types # 确保 order 与图中的顺序一致
)

# 配置并运行统计检验
# 注意：进行多重比较时，通常需要进行多重检验校正 (Multiple Testing Correction)。
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    # 多重检验校正：例如使用 Benjamini/Hochberg FDR
    comparisons_correction='fdr_bh', # <--- 关键参数！
    pvalue_thresholds=[
        [0.05, "*"],
        [0.01, "**"],
        [0.001, "***"],
        [0.0001, "****"]
    ]
)

annotator.apply_and_annotate()

plt.title('LV1 Score Pairwise Comparison by Cell Type (FDR corrected)')
plt.xlabel('Cell Type')
plt.ylabel('LV1 Score')
plt.tight_layout()
plt.show()