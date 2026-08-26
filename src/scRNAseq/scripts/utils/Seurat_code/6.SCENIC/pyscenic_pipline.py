#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Title   :   SCENIC analysis 
@File    :   pyscenic_pipline.py
@Author  :   Songqi Duan
@Contact :   songqi.duan@outlook.com
@License :   Copyright (C) 2017-2021 by Songqi Duan | 段松岐
@Created :   2021/03/07 15:39:35
@Updated :   2021/04/14 18:43:17
'''

# 导入包
import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import argparse

# 参数配置
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="count.txt")
ap.add_argument("-s", "--species", required=True,
                help="human")
                #help="human, mouse or drosophila")
ap.add_argument("-o", "--output", required=True, help="Output Folder Name")
args = vars(ap.parse_args())

INPUT_FOLDER = args["input"]
OUTPUT_FOLDER = args["output"]
SPECIES = args["species"]

if os.path.exists(OUTPUT_FOLDER) == False:
    os.makedirs(OUTPUT_FOLDER)
    print("Created New Folder")

# DATABASE_FOLDER = "/public/ref/cisTarget_databases"
DATABASE_FOLDER = "/home/chengyq/Pericyte/code/4-scenic/try-0816"
 
if SPECIES == "human":
    DATABASES_GLOB = os.path.join(
        DATABASE_FOLDER, "human", "hg38__*_tss.mc9nr.genes_vs_motifs.rankings.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(
        DATABASE_FOLDER, "human", "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
    MM_TFS_FNAME = os.path.join(
        DATABASE_FOLDER, 'resources', 'hs_hgnc_tfs.txt')
elif SPECIES == "mouse":
    DATABASES_GLOB = os.path.join(
        DATABASE_FOLDER, "mouse", "mm9-*.mc9nr.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(
        DATABASE_FOLDER, "mouse", "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
    MM_TFS_FNAME = os.path.join(DATABASE_FOLDER, 'resources', 'mm_mgi_tfs.txt')
elif SPECIES == "drosophila":
    DATABASES_GLOB = os.path.join(
        DATABASE_FOLDER, "drosophila", "dm6-*.mc9nr.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(
        DATABASE_FOLDER, "drosophila", "motifs-v8-nr.flybase-m0.001-o0.0.tbl")
    MM_TFS_FNAME = os.path.join(
        DATABASE_FOLDER, 'resources', 'allTFs_dmel.txt')
else:
    print("Incomplete parameters")

SC_EXP_FNAME = os.path.join(INPUT_FOLDER, "counts.txt")

ADJACENCIES_FNAME = os.path.join(OUTPUT_FOLDER, "adjacencies.tsv")
MODULES_FNAME = os.path.join(OUTPUT_FOLDER, "modules.p")
MOTIFS_FNAME = os.path.join(OUTPUT_FOLDER, "motifs.csv")
REGULONS_FNAME = os.path.join(OUTPUT_FOLDER, "regulons.p")
AUC_FNAME = os.path.join(OUTPUT_FOLDER, "auc.tsv")
LOOM_OUT = os.path.join(OUTPUT_FOLDER, "o.loom")

if __name__ == '__main__':
    # 读入表达矩阵，表达矩阵的格式：横坐标是基因，纵坐标是细胞
    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
    print(ex_matrix.shape)
    print(MOTIF_ANNOTATIONS_FNAME)
    print(MM_TFS_FNAME)
    # 导入转录因子
    tf_names = load_tf_names(MM_TFS_FNAME)

    # 导入数据库
    db_fnames = glob.glob(DATABASES_GLOB)

    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname))
           for fname in db_fnames]

    print(dbs)
    
    # Step 1：基于共表达情况鉴定每个TF的潜在靶点；
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
    adjacencies.to_csv(ADJACENCIES_FNAME, index=False, sep='\t')
    #adjacencies = pd.read_table(ADJACENCIES_FNAME)
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    with open(MODULES_FNAME, 'wb') as f:
        pickle.dump(modules, f)

    # Step 2: 基于motif分析选择潜在的直接结合靶点；
    # 计算所有模块的富集基序列表和相应的目标基因
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
        df.to_csv(MOTIFS_FNAME)
    
    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)
    
    # Save the enriched motifs and the discovered regulons to disk.
    with open(REGULONS_FNAME, 'wb') as f:
        pickle.dump(regulons, f)
        
    # Step 3: 分析每个细胞的regulons活性；
    regulons=[r.rename(r.name.replace('(+)', ' ('+str(len(r))+'g)')) for r in regulons]
    auc_mtx = aucell(ex_matrix, regulons, num_workers=10)
    auc_mtx.to_csv(AUC_FNAME)
