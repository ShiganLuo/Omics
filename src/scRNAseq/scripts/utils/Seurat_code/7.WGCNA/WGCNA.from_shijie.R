################################# 加载包 ################################
# rm(list = ls())
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) 
library(data.table) 
library(Seurat)
getwd()
setwd("~/brain/WGCNA/output")
### 启用WGCNA多核计算
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 
getwd()



################################# 0.输入数据准备 ################################
seurat.obj <- readRDS("~/brain/SCTed.rds")
#以样本为分析模块
seurat.obj@meta.data$sample = rownames(seurat.obj@meta.data)
seurat.obj@meta.data$sample = gsub("_\\d+$", "",seurat.obj@meta.data$sample)
saveRDS(seurat.obj, file = "~/brain/WGCNA/output/read.rds")


# # 查看表达矩阵中的细胞数量
# n_cells_counts <- ncol(seurat.obj@assays$RNA@counts)
# 
# # 查看元数据中的细胞数量
# n_cells_meta <- nrow(seurat.obj@meta.data)

# 查看当前分组
Idents(seurat.obj) <- "group"
table(Idents(seurat.obj))
#  CON  ASD   BD  MDD PTSD   SZ 
# 2719 2378 2868 1551 2513 2838


#按照样本取平均，获得表达矩阵
# 对 Seurat 对象进行对数归一化处理
seurat.obj = NormalizeData(object = seurat.obj, assay = "RNA")
expr <- AverageExpression(seurat.obj, group.by = "sample", assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
dim(expr)
# [1] 24664     50

# 更改矩阵行名为样本＋疾病，为表型数据做准备
expr2 <- as.data.frame(AverageExpression(seurat.obj, group.by = c("sample"), assays = "RNA")[["RNA"]])
name2 <- as.data.frame(colnames(expr2))
name2 <- str_split_fixed(name2$`colnames(expr2)`,"_",2)
name2 <- as.data.frame(name2)
table(name2$V1)
# ASD   BD  CON  MDD PTSD   SZ 
# 1    1    1    1    1    1 
# ASD   BD  CON  MDD PTSD   SZ 
# 16    8    4    4    4   14  


expr <- t(expr)
expr <- as.data.frame(expr)
rownames(expr) <- paste0(rownames(expr),"#",name2$V1)
expr[1:4,1:4]
#          AL627309.1 AL627309.5  LINC01409  LINC01128
# ASD#ASD 0.01300634 0.01156119 0.03528014 0.02901486
# BD#BD   0.01078420 0.03175346 0.09250293 0.06694964
# CON#CON 0.01182909 0.02654205 0.10965935 0.06104129
# MDD#MDD 0.01107854 0.04542202 0.08104490 0.07360202

#               AL627309.1 AL627309.5  LINC01409 LINC01128
# ASD_ASD1#ASD  0.03510312  0.1103397 0.24547947 0.2578503
# ASD_ASD10#ASD 0.00000000  0.0000000 0.76970140 0.6300403
# ASD_ASD11#ASD 0.00000000  0.1282051 0.07032646 0.1282051
# ASD_ASD12#ASD 0.13006622  0.0466357 0.00000000 0.1891685


### 筛选数据前7000的基因
expr <- t(expr)
keep_data <- expr[order(apply(expr, 1, mad),decreasing = T)[1:7000],]
dim(keep_data)
# [1] 7000   50
keep_data[1:4,1:4]
#           ASD_ASD1#ASD ASD_ASD10#ASD ASD_ASD11#ASD ASD_ASD12#ASD
# MALAT1    900.54503     244.30040     381.04634     372.20257
# FRMD4A     94.66663      42.20764      37.88116      63.24175
# LRMDA      53.01843      36.85925      20.21722      17.71590
# PLXDC2     61.74912      45.21459      36.47387      54.51838
keep_data <- t(keep_data)
keep_data <- as.data.frame(keep_data)
dim(keep_data)
# [1]   50 7000
keep_data[1:4,1:4]

#           PLXDC2        CALM1       KCNIP4        LSAMP
# ASD#ASD 1.091874e+14 2.950251e+01 1.163843e+04 2.329951e+31
# BD#BD   7.594853e+43 3.662962e+00 1.787717e+61 6.064693e+14
# CON#CON 1.699494e+36 3.024278e+33 4.092916e+32 2.580764e+17
# MDD#MDD 1.046094e+22 1.441165e+34 1.630571e-02 6.868630e+29

#                MALAT1   FRMD4A    LRMDA   PLXDC2
# ASD_ASD1#ASD  900.5450 94.66663 53.01843 61.74912
# ASD_ASD10#ASD 244.3004 42.20764 36.85925 45.21459
# ASD_ASD11#ASD 381.0463 37.88116 20.21722 36.47387
# ASD_ASD12#ASD 372.2026 63.24175 17.71590 54.51838

### 创建datTraits，包含分组、表型等信息
df <- data.frame(rownames(keep_data))
a <- as.data.frame(str_split_fixed(df$rownames.keep_data.,"#",2))
df <- cbind(df,a)
colnames(df) <- c("AverageExpression_group","Sample","group")

datTraits <- data.frame(row.names = rownames(keep_data),
                        group = df$group)
table(datTraits$group)
# ASD   BD  CON  MDD PTSD   SZ 
#  1    1    1    1    1    1 
# ASD   BD  CON  MDD PTSD   SZ 
# 16    8    4    4    4   14 

### 给分组加上编号
grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
head(datTraits)
#            group groupNo
# ASD#ASD     ASD       1
# BD#BD        BD       2
# CON#CON     CON       3
# MDD#MDD     MDD       4
# PTSD#PTSD  PTSD       5
# SZ#SZ        SZ       6
#                 group groupNo
# ASD_ASD1#ASD    ASD       1
# ASD_ASD10#ASD   ASD       1
# ASD_ASD11#ASD   ASD       1
# ASD_ASD12#ASD   ASD       1
# ASD_ASD13#ASD   ASD       1
# ASD_ASD14#ASD   ASD       1

table(datTraits$group)
table(datTraits$groupNo)

#使用输入文件
datExpr0 <- as.data.frame(keep_data)




############################## 1.判断数据质量 ################################
### 判断数据质量--缺失值
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){ #如果存在异常样本或基因
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) # 异常的基因
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) # 异常的样本
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:从数据中删除违规基因和样本
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

### 绘制样品的系统聚类树
if(T){
  #针对样本做聚类树
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  pdf("step1_Sample dendrogram.pdf",width = 8,height = 6)
  p <- plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
            cex.axis = 1, cex.main = 1,cex.lab=1)
  print(p)
  dev.off()
  # ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)),
                                  colors = rainbow(length(table(datTraits$group))),
                                  signed = FALSE)
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6)
  p2 <- plotDendroAndColors(sampleTree, sample_colors,
                            groupLabels = "trait",
                            cex.dendroLabels = 0.8,
                            marAll = c(1, 4, 3, 1),
                            cex.rowText = 0.01,
                            main = "Sample dendrogram and trait"
  )
  print(p2)
  dev.off()
  # # Plot a line to show the cut
  # abline(h = 23500, col = "red") #根据实际情况而定
  # dev.off()
}
###____________________________________________________________________________________________________________________
##若存在显著离群点；剔除掉
if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 1500, minSize = 10) # cutHeight根据实际情况而定
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}

### 判断数据质量 : PCA进行分组查看
group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, #标签不重叠
                    col.ind = group_list, # 分组上色
                    axes.linetype=NA,  # remove axeslines
                    mean.point=F#去除分组中心点
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) #坐标轴的纵横比
pca
ggsave(pca,filename= "step1_Sample PCA analysis.pdf", width = 8, height = 8)

##保存数据
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")



############################### 2.挑选最佳阈值power ###################################
rm(list = ls())  
load("step1_input.Rdata")
R.sq_cutoff = 0.8  #设置R^2 cut-off，默认为0.85
if(T){
  # Call the network topology analysis function
  #设置power参数选择范围
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  pdf("step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: 寻找拐点，确认最终power取值
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate  #查看估计的最佳power
# [1] 4
power = sft$powerEstimate
power 
# [1] 4

# 若无向网络在power小于15或有向网络powerdx小于30内，没有一个power值使
# 无标度网络图谱结构R^2达到0.8且平均连接度在100以下，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if(is.na(power)){
  # 官方推荐 "signed" 或 "signed hybrid"
  # 为与原文档一致，故未修改
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

save(sft, power, file = "step2_power_value.Rdata")


##################### 3.一步法构建加权共表达网络，识别基因模块 ####################
rm(list = ls())  
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr), #默认5000
    corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 70,    ##越大模块越少 30
    mergeCutHeight = 0.25, ##越大模块越少 0.25
    numericLabels = TRUE, 
    saveTOMs = T,
    verbose = 3
  )
  table(net$colors) 
  # 0    1    2    3    4    5    6    7    8    9   10   11   12   13 
  # 2389 1000  662  499  439  359  340  280  278  225  158  139  117  115 
  
  # power: 上一步计算的软阈值
  # maxBlockSize:计算机能处理的最大模块的基因数量(默认5000),16G内存可以处理2万个，
  # 计算资源允许的情况下最好放在一个block里面。
  # corType：计算相关性的方法；可选pearson(默认)，bicor。后者更能考虑离群点的影响。
  # networkType:计算邻接矩阵时，是否考虑正负相关性；默认为"unsigned",可选"signed", "signed hybrid"
  # TOMType：计算TOM矩阵时，是否考虑正负相关性；默认为"signed",可选"unsigned"。但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义。
  # numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  # saveTOMs：最耗费时间的计算，可存储起来供后续使用，
  # mergeCutHeight: 合并模块的阈值，越大模块越少,一般为0.25
  # minModuleSize: 每个模块里最少放多少个基因，设定越大模块越少
  # 输出结果根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
  # **0 (grey)**表示**未**分入任何模块的基因。
}


## 模块可视化，层级聚类树展示各个模块
if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "step3_genes_modules.Rdata")


################################# 4.关联基因模块与表型 #####################################
rm(list = ls())
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
datTraits$group
# datTraits$group <- ifelse(grepl("Normal#$", rownames(datTraits)), "Normal", "COVID-19")
## 模块与表型的相关性热图
if(T){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group) #get the group
  
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("step4_Module-trait-relationship_heatmap2.pdf",
      width = 1.5*length(colnames(design)),
      height = 0.4*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) #留白：下、左、上、右
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 # colors = blueWhiteRed(50),
                 colors = colorRampPalette(c("#105eb7", "white","#d7131a"))(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = "step4_design.Rdata")
}
#条形图
if(T){
  y = as.data.frame(design[ ,"SZ"]);
  GS1=as.numeric(cor(y,datExpr,use="p"))
  GeneSignificance=abs(GS1)
  # Next module significance is defined as average gene significance.
  ModuleSignificance=tapply(GeneSignificance,moduleColors , mean , na.rm=T)
  pdf("step4-Module-single_trait-relationship_barplot_SZ.pdf",width = 8,height = 5)
  plotModuleSignificance(GeneSignificance,moduleColors)
  dev.off()
}

### 模块与表型的相关性boxplot图
if(T){
  mes_group <- merge(MEs,datTraits,by="row.names")
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = c("#CC0066","#99FFFF", "#99FF33","#F39B7FCC","#FFFF99","#CCCCFF"),
              #palette = "jco",
              #add="jitter",
              legend = "") +stat_compare_means()
  }  
  # 批量画boxplot
  colorNames <- names(MEs)
  pdf("step4_Module-trait-relationship_boxplot.pdf", width =16,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=4)) #排布为每行4个
  dev.off()
}


### 基因与模块、表型的相关性散点图
#所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因算出相关系数，
#如果跟性状显著相关的基因也跟某个模块显著相关，那么这些基因可能就非常重要。

# 选择离散性状的表型
levels(datTraits$group)
# [1] "ASD"  "BD"   "CON"  "MDD"  "PTSD" "SZ" 
choose_group <- "CON"

if(T){
  modNames <- substring(names(MEs), 3)
  
  ### 计算模块与基因的相关性矩阵
  ## Module Membership: 模块内基因表达与模块特征值的相关性
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  
  ###  计算性状与基因的相关性矩阵
  ## Gene significance，GS：比较样本某个基因与对应表型的相关性
  ## 连续型性状
  # trait <- datTraits$groupNo
  ## 非连续型性状，需转为0-1矩阵, 已存于design中
  trait <- as.data.frame(design[,choose_group])
  geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  ### 可视化基因与模块、表型的相关性.
  #selectModule<-c("blue","green","purple","grey")  ##可以选择自己想要的模块
  selectModule <- modNames  ## 全部模块批量作图
  pdf("step4_gene-Module-trait-significance_CON.pdf",width=7, height=1.5*ncol(MEs))
  par(mfrow=c(ceiling(length(selectModule)/2),2)) #批量作图开始
  for(module in selectModule){
    column <- match(module,selectModule)
    print(module)
    moduleGenes <- moduleColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for trait",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}


#########################  5. WGCNA可视化：TOMplot  Eigengene-adjacency-heatmap ##################################
rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

if(T){
  TOM=TOMsimilarityFromExpr(datExpr,power=power)
  dissTOM=1-TOM
  ## draw all genes
  if(T){
    geneTree = net$dendrograms[[1]]
    plotTOM = dissTOM^7
    diag(plotTOM)=NA
    png("step5_TOMplot_Network-heatmap.png",width = 800, height=600)
    TOMplot(plotTOM,geneTree,moduleColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot")
    dev.off()
  }
  ### draw selected genes to save time...just for test...
  if(T){
    nSelect =0.1*nGenes
    set.seed(123)
    select=sample(nGenes,size = nSelect)
    selectTOM = dissTOM[select,select]
    selectTree = hclust(as.dist(selectTOM),method = "average")
    selectColors = moduleColors[select]
    plotDiss=selectTOM^7
    diag(plotDiss)=NA
    pdf("step5_select_TOMplot_Network-heatmap.pdf",width=8, height=6)
    TOMplot(plotDiss,selectTree,selectColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot of selected gene")
    dev.off()
  }
}


### 模块相关性展示 Eigengene-adjacency-heatmap
if(T){
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # 若添加表型数据
  if(T){
    ## 连续型性状
    # MET = orderMEs(cbind(MEs,datTraits$groupNo))
    # 非连续型性状，需将是否属于这个表型进行0,1数值化，已存于design中
    # design
    MDD = as.data.frame(design[,1])
    names(MDD) = "MDD"
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, MDD))
  }
  pdf("step5_module_cor_Eigengene-dendrogram_MDD.pdf",width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="",
                        marDendro = c(0,4,1,4),  # 留白：下右上左
                        marHeatmap = c(5,5,1,2), # 留白：下右上左
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}

#################### 6. 选择感兴趣基因模块进行GO分析 ####################
rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

### 条件设置
OrgDb = "org.Hs.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"
table(moduleColors)
# moduleColors
# black        blue       brown       green greenyellow        grey     magenta 
# 280         662         499         359         139        2389         225 
# pink      purple         red      salmon         tan   turquoise      yellow 
# 278         158         340         115         117        1000         439 

choose_module <- c("black","blue","brown","green","greenyellow", "grey","magenta"
                   ,"pink","purple","red","salmon","tan","turquoise","yellow") 

if(T){
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  
  gene_module <- data.frame(gene=colnames(datExpr),
                            module=moduleColors)
  write.csv(gene_module,file = "step6_gene_moduleColors.csv",row.names = F, quote = F)
  tmp <- bitr(gene_module$gene,fromType = genetype,  # "SYMBOL"   "ENSEMBL"
              toType = "ENTREZID",
              OrgDb = OrgDb )
  gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  ###run go analysis
  formula_res <- compareCluster(
    ENTREZID~module,
    data = choose_gene_module_entrz,
    fun = "enrichGO",
    OrgDb = OrgDb,
    ont = "BP",  #One of "BP", "MF", and "CC"  or "ALL"
    pAdjustMethod = "BH",   #指定多重假设检验矫正的方法,一般选择 "BH" 或 "fdr"，BH较严格，fdr较温和（计算的q小些）
    pvalueCutoff = 0.05
    # qvalueCutoff = 0.25
    # minGSSize &maxGSSize：是富集的最小/大的基因集的大小（基因数目）
  )
  
  ###精简GO富集的结果,去冗余
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res,
    cutoff=0.5,
    by="p.adjust",
    select_fun=min
  )
  save(gene_module, formula_res, lineage1_ego, file="step6_module_GO_term.Rdata")
  write.csv(lineage1_ego@compareClusterResult,
            file="step6_module_GO_term.csv")
  ### 绘制dotplot图
  dotp <- dotplot(lineage1_ego,
                  showCategory=10,
                  includeAll = TRUE, #将有overlap的结果也展示出来
                  label_format=90)
  ggsave(dotp,filename= "step6_module_GO_term.pdf", #device = cairo_pdf,
         width = 16,
         height = 15)
}

############################### 7.感兴趣基因模块绘制热图 ######################################
rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step3_genes_modules.Rdata")
table(moduleColors)
# moduleColors
# black        blue       brown       green greenyellow        grey     magenta 
# 280         662         499         359         139        2389         225 
# pink      purple         red      salmon         tan   turquoise      yellow 
# 278         158         340         115         117        1000         439 

module = "blue"
### 感兴趣模块画热图
if(T){
  dat=datExpr[,moduleColors==module]
  library(pheatmap)
  n=t(scale(dat)) #对基因做scale，并转置表达矩阵为行为基因、列为样本形式
  # n[n>2]=2
  # n[n< -2]= -2
  # n[1:4,1:4]
  
  group_list=datTraits$group
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  pheatmap::pheatmap(n,
                     fontsize = 8,
                     show_colnames =T,
                     show_rownames = F,
                     cluster_cols = T,
                     annotation_col =ac,
                     width = 8,
                     height = 6,
                     angle_col=45,
                     main = paste0("module_",module,"-gene heatmap"),
                     filename = paste0("step7_module_",module,"_Gene-heatmap.pdf"))
  
}

################### 8.感兴趣模块基因导出 VisANT or cytoscape ######################
rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")
module = "red"  
if(T){
  ### 提取感兴趣模块基因名
  gene <- colnames(datExpr)
  inModule <- moduleColors==module
  modgene <- gene[inModule]
  # write.table(modgene,paste0("step8_",module,"_modgene.csv"))
  
  ### 模块对应的基因关系矩阵
  TOM <- TOMsimilarityFromExpr(datExpr,power=power)
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modgene,modgene)
  
  ### 筛选连接度最大的top100基因（核心基因/hub基因）
  nTop = 100#30
  IMConn = softConnectivity(datExpr[, modgene]) #计算连接度
  top = (rank(-IMConn) <= nTop) #选取连接度最大的top100
  filter_modTOM <- modTOM[top, top]
  
  # for visANT
  vis <- exportNetworkToVisANT(filter_modTOM,
                               file = paste("step8_visANTinput-",module,".txt",sep = ""),
                               weighted = T,threshold = 0)
  # for cytoscape
  cyt <- exportNetworkToCytoscape(filter_modTOM,
                                  edgeFile = paste("step8_CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("step8_CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                  weighted = TRUE,
                                  threshold = 0.15,  #weighted权重筛选阈值，可调整
                                  nodeNames = modgene[top],
                                  nodeAttr = moduleColors[inModule][top])
}


#######WGCNA种提取模块gene进行打分-----
modulegene <- read.csv("~/brain/WGCNA/output/step6_gene_moduleColors.csv")
modulegene <- subset(modulegene,subset = module %in% "red") %>% pull(gene)
library(ggplot2)
library(ggpubr)
####模块基因进行功能富集-----
library(DOSE)
library(enrichplot)
library(clusterProfiler)
modulegene <- read.csv("~/brain/WGCNA/output/step6_gene_moduleColors.csv")
red <- subset(modulegene,subset = module %in% "red") %>% pull(gene)

library(org.Hs.eg.db)
bp <-
  enrichGO(
    red,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term <- bp@result
saveRDS(term,"~/brain/WGCNA/output/term_red.rds")
library(ggplot2)
# Cluster2
#term <- readRDS("~/brain/WGCNA/output/term.rds")
go <- c("dendrite morphogenesis","dendrite development","dendritic spine morphogenesis",
        "dendritic spine development","regulation of postsynapse organization","dendritic spine organization",
        "vesicle organization","neuron projection organization","postsynapse organization","dendritic transport","neurotrophin signaling pathway")

##brown
# "myelination","ensheathment of neurons","axon ensheathment","glial cell differentiation",
# "regulation of myelination","sphingolipid metabolic process","gliogenesis","regulation of ubiquitin-dependent protein catabolic process",
# "sphingolipid biosynthetic process","myelin assembly","fatty acid elongation","substantia nigra development",
# "radial glial cell differentiation","neuron projection arborization"


# black:
# "glial cell differentiation","cholesterol biosynthetic process","oligodendrocyte differentiation",
# "gliogenesis","regulation of glial cell differentiation","axon development",
# "central nervous system myelination","axon ensheathment in central nervous system","ensheathment of neurons",
# "axon ensheathment","oligodendrocyte development","glial cell development","regulation of gliogenesis",
# "regulation of neuron projection development","myelination","regulation of nervous system development","astrocyte differentiation","axonogenesis"


term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
# pdf("~/20230410_heart-SMC/result/resultdim30_newanno/Figure6-monocle-new/monocle2/1e-9/Figure/分支基因/figure/1.pdf",height = 2,width = 4)
p <- ggplot(data = term,
            aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "red",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=3.5,hjust =0,color="white")+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
print(p)
dev.off()
