# Description

## DESeq.r

for DESeq.r, gsea.r you can Rscript script.r --help for quickly know how to use it.

20250324:

    - 修复：基因背景下的TE差异分析，绘制火山图和热图时需要把基因和TE单独提出来进行绘制

20250326:

    - 修复TEcoutnmode
    - DESeq2.r支持根据group.csv调整DESeq2condition
    - 修复plot_PCANorm写死样本数目和样本组名称的错误
    - plot_heatmap()函数添加判断只有上调、下调或失调基因数量大于2才进行热图绘制，防止程序中断

## geas.r

for gsea.r, gsea.r you can Rscript script.r --help for quickly know how to use it.
to do: fulfill plot_enrichment for beautiful plot (now just use pygsea to achieve it)

## go-kegg.r

temporary for mouse

for go-kegg.r you can Rscript script.r --help for quickly know how to use it.

20250326:
    - 支持命令行方式传递参数
    - 合并kegg.r和go.r
