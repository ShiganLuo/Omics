import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import colorsys
import numpy as np


def genotype(text:str):
# match genotype
# print(text)
    #有的基因型以|分割，且在相位中也可能展示基因型；所以第一个\t和:包括内容才为基因型
    pattern1 = re.compile(r'\t(\d.\d):')
    pattern2 = re.compile(r'(\d*\t)(\d*\t)(.\t)([a-zA-Z]*\t)([a-zA-Z]*\t)')
    genotype = pattern1.findall(text)
    match = pattern2.match(text)
    print(match.groups())
    w1 = match.groups()[3].replace("\t","")
    w2 = match.groups()[4].replace("\t","")
    print(w1)
    print(w2)
    # print(genotype)
    # print(len(genotype))
    # plot
    #w1:0；w2:1
    g1 = genotype.count('0/0') + genotype.count('0|0')
    g2 = genotype.count('0/1') + genotype.count('1/0') + genotype.count('0|1') + genotype.count('1|0')
    g3 = genotype.count('1/1') + genotype.count('1|1')
    print(f'{w1}{w1}:{g1}')
    print(f'{w1}{w2}:{g2}')
    print(f'{w2}{w2}:{g3}')
    sta = {f'{w1}/{w1}':g1,f'{w1}/{w2}':g2,f'{w2}/{w2}':g3}
    return sta

def genetotype_pie(texts:list,outdir:str):
    ####饼图
    fig,(ax1,ax2,ax3) = plt.subplots(1, 3, sharex=True,figsize=(10, 6))
    Domestication = genotype(texts[1])
    print(Domestication)
    Dvalue = list(Domestication.values())
    Dlabel = list(Domestication.keys())
    Splendens = genotype(texts[3])
    print(Splendens)
    Svalue = list(Splendens.values())
    Sibling = genotype(texts[5])
    Sivalue = list(Sibling.values())
    print(Sibling)
    ax1.pie(Dvalue,colors=['red','green','blue'])
    ax2.pie(Svalue,colors=['red','green','blue'])
    ax3.pie(Sivalue,colors=['red','green','blue'])
    fig.legend(Dlabel, loc="upper right")
    #标记需手动修改
    ax1.text(-0.4,1.2,"Domestication")
    ax1.text(-0.1,1.4,"479")
    ax2.text(-0.35,1.2,"Splendens")
    ax2.text(-0.1,1.4,"27")
    ax3.text(-0.35,1.1,"Relatives")
    ax3.text(-0.1,1.4,"93")
    fig.text(0.45,0.05,"Chr21:10201624")
    fig.savefig(f'{outdir}/PLA2G1B_Chr21_10201624.jpg')

def vary_hue(color, factor):
    """在父类颜色基础上调整色相，alpha固定1"""
    r,g,b = color[:3]
    h,l,s = colorsys.rgb_to_hls(r,g,b)
    h = (h + factor) % 1.0
    r_new,g_new,b_new = colorsys.hls_to_rgb(h,l,s)
    return (r_new,g_new,b_new,1.0)

def resolve_label_positions(label_positions, pad=0.02):
    """避免标签重叠"""
    left,right=[],[]
    for ang, tx, ty, text in label_positions:
        (right if np.cos(np.deg2rad(ang))>=0 else left).append([ang,tx,ty,text])

    left.sort(key=lambda x:-x[2])
    right.sort(key=lambda x:-x[2])

    def adjust(side):
        if not side: return side
        for i in range(1,len(side)):
            prev_y = side[i-1][2]
            curr_y = side[i][2]
            if curr_y > prev_y - pad: side[i][2] = prev_y - pad
        for i in range(len(side)-2,-1,-1):
            next_y = side[i+1][2]
            curr_y = side[i][2]
            if curr_y < next_y + pad: side[i][2] = next_y + pad
        for i in range(len(side)):
            side[i][2] = min(max(side[i][2],0.02),0.98)
        return side

    return adjust(left) + adjust(right)

def draw_arc_text_quadrant(ax, center, radius, text, theta1, theta2, **kwargs):
    """
    内层文字沿扇区切线旋转，按照扇区象限区分方向
    """
    theta_mid = (theta1 + theta2) / 2
    x = center[0] + radius * np.cos(np.deg2rad(theta_mid))
    y = center[1] + radius * np.sin(np.deg2rad(theta_mid))

    # 根据象限调整旋转角度
    if 0 <= theta_mid <= 90:  # 右上
        rotation = theta_mid - 90
    elif 90 < theta_mid <= 180:  # 左上
        rotation = theta_mid - 90
    elif 180 < theta_mid <= 270:  # 左下
        rotation = theta_mid + 90
    else:  # 右下 270~360
        rotation = theta_mid - 90

    ax.text(x, y, text, rotation=rotation, ha='center', va='center', **kwargs)



def nested_donut_discrete(df: pd.DataFrame, columns: list, outplot:str, figsize=(8,8), startangle=90):
    """内层标签沿弧完整显示且方向自然，外层标签避让，父子颜色离散"""
    # 按类别数排序（内层少->外层多）
    columns_sorted = sorted(columns, key=lambda c: df[c].nunique())
    layer_values, layer_labels, layer_colors = [], [], []

    # 内层父类颜色（离散、饱和）
    col0 = columns_sorted[0]
    group0 = df.groupby(col0, sort=False).size()
    layer_values.append(group0.values.tolist())
    layer_labels.append(group0.index.tolist())
    n0 = len(group0)
    base_colors = sns.color_palette("tab20", n0)
    layer_colors.append([(*c,1.0) for c in base_colors])

    # 外层子类颜色
    for layer_idx in range(1, len(columns_sorted)):
        col = columns_sorted[layer_idx]
        parent_col = columns_sorted[layer_idx-1]
        group = df.groupby([parent_col, col], sort=False).size()
        outer_parents = group.index.get_level_values(0).tolist()
        outer_labels = group.index.get_level_values(1).tolist()
        outer_values = group.values.tolist()
        colors = []
        parent_index = {p:i for i, p in enumerate(layer_labels[layer_idx-1])}
        parent_seen = {}
        for p in outer_parents:
            idx = parent_index[p]
            base = layer_colors[layer_idx-1][idx]
            count = outer_parents.count(p)
            k = parent_seen.get(p,0)
            factor = -0.2 + 0.4*k/max(count-1,1)
            c = vary_hue(base, factor)
            colors.append(c)
            parent_seen[p] = k+1
        layer_labels.append(outer_labels)
        layer_values.append(outer_values)
        layer_colors.append(colors)

    # 绘图
    fig = plt.figure(figsize=figsize)
    ax_main = fig.add_axes([0.15,0.15,0.7,0.7])
    ax_main.axis('equal')
    ax_label = fig.add_axes([0.02,0.02,0.96,0.96])
    ax_label.axis('off')

    radius = 1.0
    width = 0.28
    outer_info = []

    for layer_idx, (values, labels, colors) in enumerate(zip(layer_values, layer_labels, layer_colors)):
        wedges,_ = ax_main.pie(values, radius=radius, colors=colors, labels=None,
                               startangle=startangle, wedgeprops=dict(width=width, edgecolor='white'))
        # 内层标签沿圆弧且方向自然
        if layer_idx == 0:
            for w, lab in zip(wedges, labels):
                draw_arc_text_quadrant(ax_main, (0,0), radius-width/2, lab, w.theta1, w.theta2,
                                    fontsize=10, color='black')


        # 外层记录用于标签
        if layer_idx == len(columns_sorted)-1:
            for w, lab in zip(wedges, labels):
                ang = (w.theta2 + w.theta1)/2
                outer_info.append((ang, lab, w))

        radius += width*0.98  # 紧贴

    # 外层标签
    raw_positions = []
    for ang, text, w in outer_info:
        x = np.cos(np.deg2rad(ang))
        y = np.sin(np.deg2rad(ang))
        px = 0.5 + x*0.35*1.05
        py = 0.5 + y*0.35*1.05
        tx = px + 0.10*x
        ty = py + 0.10*y
        tx = min(max(tx,0.02),0.98)
        ty = min(max(ty,0.02),0.98)
        raw_positions.append([ang, tx, ty, text])

    final_positions = resolve_label_positions(raw_positions, pad=0.02)

    for ang, tx, ty, text in final_positions:
        x = np.cos(np.deg2rad(ang))
        y = np.sin(np.deg2rad(ang))
        px = 0.5 + x*0.35*1.05
        py = 0.5 + y*0.35*1.05
        ax_label.text(tx, ty, text, ha='left' if x>=0 else 'right', va='center',
                      fontsize=10, transform=ax_label.transAxes)
        ax_label.plot([px, tx],[py, ty], color='gray', lw=0.8, transform=ax_label.transAxes)

    plt.savefig(outplot, dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    # reda file 
    infile = "/home/lsg/GWAS/Betta_splendens240901/gene_Investication/output/snp/dip2a.txt"
    outdir = "/home/lsg/GWAS/Betta_splendens240901/gene_Investication/output/plot/pie"
    fo = open(infile,'r')
    texts = fo.readlines()
    fo.close()
    genetotype_pie(texts,"a.png")

    # two degree nested_donut
    df = pd.read_csv("/home/luosg/Data/genomeStability/data/target_fq_deNormal.tsv",sep="\t")
    nested_donut_discrete(df, ["Status","Cell_line"],"/home/luosg/Data/genomeStability/output/result/metadata/nested_donut_discrete.png")

