import pandas as pd
from collections import Counter

df1 = pd.read_csv("/disk5/luosg/Totipotent20251031/data/geneset/8C1.tsv",header=None)
df2 = pd.read_csv("/disk5/luosg/Totipotent20251031/data/geneset/8C2.tsv",header=None)
df3 = pd.read_csv("/disk5/luosg/Totipotent20251031/data/geneset/8C3.tsv",header=None)

# 取三列的第一列
g1 = df1.iloc[:,0].astype(str)
g2 = df2.iloc[:,0].astype(str)
g3 = df3.iloc[:,0].astype(str)

# 计数
cnt = Counter(list(g1) + list(g2) + list(g3))

# 取出现 >= 2 次的基因
common = sorted([g for g, n in cnt.items() if n >= 2])
geneset_name = "act_8cell_gene"
description = "Homo sapiens"

gmt_line = geneset_name + "\t" + description + "\t" + "\t".join(common)

with open("/disk5/luosg/Totipotent20251031/data/geneset/8C_human.gmt", "w") as f:
    f.write(gmt_line + "\n")