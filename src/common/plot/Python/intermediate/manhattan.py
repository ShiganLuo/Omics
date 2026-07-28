# %%
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
import re
from pathlib import Path
import sys
colorList=['#11604B','#13234A']
#################################-------------------------------------函数----------------------------------###############
# mpl.rcParams["font.sans-serif"] = ["SimHei"]
# # 设置正常显示符号
# mpl.rcParams["axes.unicode_minus"] = False

# print(mpl.matplotlib_fname())

my_font = FontProperties(fname=r"/home/lsg/anaconda3/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/SimHei.ttf",size=36)
# %%
#绘制曼哈顿图主体函数
def manhattan(data,chr,pos,value,ax,palette='Set2',log=True,**arg):
    r'''
    汇总绘制 Manhattan 图

    Parameters
    ----------
    data : pd.DataFrame
        数据框，包含染色体、碱基位点、值，要求已经按照染色体排序
    chr : str
        染色体
    pos : str
        碱基在染色体上的位置
    '''
    chr_list = ["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9",
        "Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21"]
  
    data = data[[chr,pos,value]].copy()
    # 计算 -log10，log参数控制是否进行负对数转换
    data.loc[:,'-logp'] = -np.log10(data[value]) if log else data[value]
  
    tmp = data.groupby(chr,sort=False)[pos].max()#每组染色体序列最大值
    tmp = tmp.cumsum()-tmp#计算这些最大值的累积和，即从第一个染色体到当前染色体的所有最大值之和。
    data[pos] = data[pos] + tmp[data[chr]].values#更新染色体绘制位置
    data = data.sort_values(by=pos)#避免存在染色体位置不是按顺序排列的情况，以防杂线出现

    mark_tmp = data.loc[data['-logp']>=-np.log10(1/19878711),:]
    for chrname in data[chr].unique():
        mark_tmp_chr = mark_tmp[mark_tmp[chr]==chrname]
        tmp_chr = data[data[chr]==chrname]
        color = palette[chr_list.index(chrname)]
        ax.plot(tmp_chr[pos],tmp_chr['-logp'],color=color)#log为false，logp列和xpclr_norm列一致
        # ax.scatter(mark_tmp_chr[pos],mark_tmp_chr['-logp'],color=color,s=250)#突出显示显著性点
    # ax.scatter(mark_tmp[pos],mark_tmp['-logp'],c='red',s=500)

    chrom_df=data.groupby(chr,sort=False)[pos].median()#计算染色体位置中位数
    ax.set_xticks(chrom_df)
    ax.set_xticklabels(chrom_df.index)
    ax.set_ylabel(f'-log({value})' if log else value,fontsize='x-large')
    ax.set_xlim(0-data[pos].max()*0.01,data[pos].max()*1.01)

#################################-------------------------------------执行----------------------------------###############
# %%
#读取文件
out = Path('plot')
f = Path(sys.argv[1])

# title = f.name[:-10]
title = f.name.split('_')[2].split('.')[0]
title = title + "_fst"
print(title)

df = pd.read_table(f)
#筛选chr列，过滤掉一些碎片序列
chrs = ["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9",
        "Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21"]
df = df[df.CHROM.isin(chrs)]
chrs_map = {chrs[i]:i+1 for i in range(len(chrs))}#创建映射染色体字典 

df['chrnumid'] = df.CHROM.map(chrs_map)#为数据框对应染色体编号添加chrnumid列 
df = df[df['MEAN_FST'] > 0]

df = df.sort_values('chrnumid')

#创建图像与定义字体大小
mpl.rcParams['font.size'] = 25#统一字体大小
fig = plt.figure(figsize=(45,10),dpi=300)
ax:plt.Axes = fig.subplots()
# y 轴标签字体大小
ax.tick_params(axis='y',labelsize='large')

#绘图
manhattan(df,'CHROM','BIN_START','MEAN_FST',ax,s=500,log = False,palette=sns.color_palette(colorList,n_colors=len(chrs)))

#细节设置
#控制y轴范围与显著线条
ax.axhline(0.36,ls='--',c='black')
ax.axhline(0.4,ls='-',c='black')
ax.set_title(title.replace('_',' '),fontproperties=my_font)
# ax.set_ylim(0,ax.get_ylim()[1])#控制负数是否显示
# ax.set_ylabel('$-log(xpclr_{norm})$',fontsize='x-large')
_ = ax.set_xticklabels(ax.get_xticklabels(),rotation=45)

ax.tick_params(axis='y',labelsize='x-large')

ax.set_xlabel('')
# break
fig.savefig(out / f'{title}.jpg',bbox_inches='tight')


