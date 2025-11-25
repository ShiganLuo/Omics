import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def multilayer_donut(df, level_cols, count_col="Count"):
    """
    df: 包含层级分类的 DataFrame
    level_cols: 按顺序从父类到子类，如 ["Level1", "Level2", "Level3"]
    """

    # 每层分类的唯一组合
    layers = []
    for i, col in enumerate(level_cols):
        group = df.groupby(level_cols[:i+1])[count_col].sum().reset_index()
        layers.append(group)

    # 基础画布
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis("equal")

    # 颜色方案：为最内层父类分配颜色，然后按子类分支派生
    base_colors = plt.cm.tab20(np.linspace(0, 1, layers[0].shape[0]))
    color_map = {}

    # 给 Level1 生成颜色
    for i, key in enumerate(layers[0][level_cols[0]]):
        color_map[(key,)] = base_colors[i]

    # 给其他层生成颜色（基于父类颜色加深/变浅）
    for level in range(1, len(level_cols)):
        parent_cols = level_cols[:level]   # 父类列
        curr_cols = level_cols[:level+1]   # 当前层列

        for _, row in layers[level].iterrows():
            parent_key = tuple(row[parent_cols])
            curr_key = tuple(row[curr_cols])

            # 父类颜色
            parent_color = color_map[parent_key]

            # 生成子类颜色（明度变化）
            rgb = mcolors.to_rgb(parent_color)
            factor = 0.6 + 0.4*np.random.rand()  # 随机在父色基础变亮或变暗
            child_color = tuple([min(1, c*factor) for c in rgb])

            color_map[curr_key] = child_color

    # 绘图，从内到外
    radius = 0.3
    width = 0.3

    for level in range(len(level_cols)):
        layer = layers[level]
        labels = []
        sizes = []
        colors = []

        for _, row in layer.iterrows():
            key = tuple(row[level_cols[:level+1]])
            labels.append(key[-1])
            sizes.append(row[count_col])
            colors.append(color_map[key])

        wedges, _ = ax.pie(
            sizes,
            radius=radius + width*level,
            colors=colors,
            startangle=90,
            counterclock=False,
            wedgeprops=dict(width=width, edgecolor='white')
        )

        # 只给最外层加标签
        if level == len(level_cols) - 1:
            for w, lab in zip(wedges, labels):
                ang = (w.theta2 + w.theta1) / 2.
                x = (radius + width*(level+0.5)) * np.cos(np.deg2rad(ang))
                y = (radius + width*(level+0.5)) * np.sin(np.deg2rad(ang))
                ax.text(x, y, lab, ha="center", va="center", fontsize=10)

    fig.savefig("a.png")




if __name__ == "__main__":
    df = pd.read_csv("/home/luosg/Data/genomeStability/data/target_fq_deNormal.tsv",sep="\t")
    multilayer_donut(df,["Status","Cell_line"])
