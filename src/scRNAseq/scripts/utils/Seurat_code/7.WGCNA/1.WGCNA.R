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
setwd("~/精神疾病/data/388/Seurat/data/WGCNA")
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 
getwd()



seurat.obj@meta.data$sample = rownames(seurat.obj@meta.data)
seurat.obj@meta.data$sample = gsub("_\\d+$", "",seurat.obj@meta.data$sample)
saveRDS(seurat.obj, file = "~/精神疾病/data/388/Seurat/data/WGCNA/WGCNA_sample.rds")

seurat.obj <- readRDS("~/精神疾病/data/388/Seurat/data/WGCNA/WGCNA_sample.rds")




Idents(seurat.obj) <- "group"
table(Idents(seurat.obj))


seurat.obj = NormalizeData(object = seurat.obj, assay = "RNA")
expr <- AverageExpression(seurat.obj, group.by = "sample", assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
dim(expr)



expr2 <- as.data.frame(AverageExpression(seurat.obj, group.by = c("sample"), assays = "RNA")[["RNA"]])
name2 <- as.data.frame(colnames(expr2))
name2 <- str_split_fixed(name2$`colnames(expr2)`,"_",2)
name2 <- as.data.frame(name2)
table(name2$V1)


expr <- t(expr)
expr <- as.data.frame(expr)
rownames(expr) <- paste0(rownames(expr),"#",name2$V1)
expr[1:4,1:4]



expr <- t(expr)
keep_data <- expr[order(apply(expr, 1, mad),decreasing = T)[1:7000],]
dim(keep_data)
keep_data[1:4,1:4]
keep_data <- t(keep_data)
keep_data <- as.data.frame(keep_data)
dim(keep_data)
keep_data[1:4,1:4]



df <- data.frame(rownames(keep_data))
a <- as.data.frame(str_split_fixed(df$rownames.keep_data.,"#",2))
df <- cbind(df,a)
colnames(df) <- c("AverageExpression_group","Sample","group")

datTraits <- data.frame(row.names = rownames(keep_data),
                        group = df$group)
table(datTraits$group)

grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
head(datTraits)

table(datTraits$group)
table(datTraits$groupNo)

datExpr0 <- as.data.frame(keep_data)




gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

if(T){
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  pdf("step1_Sample dendrogram.pdf",width = 8,height = 6)
  p <- plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
            cex.axis = 1, cex.main = 1,cex.lab=1)
  print(p)
  dev.off()
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)),
                                  colors = rainbow(length(table(datTraits$group))),
                                  signed = FALSE)
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
}

if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 1500, minSize = 10)
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}

group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE,
                    col.ind = group_list,
                    axes.linetype=NA,
                    mean.point=F#去除分组中心点
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1)
pca
ggsave(pca,filename= "step1_Sample PCA analysis.pdf", width = 8, height = 8)

datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")



rm(list = ls())  
load("step1_input.Rdata")
R.sq_cutoff = 0.8
if(T){
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  pdf("step2_power-value.pdf",width = 16,height = 12)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=R.sq_cutoff ,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate
power = sft$powerEstimate
power 

if(is.na(power)){
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


rm(list = ls())  
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 45,
    mergeCutHeight = 0.25,
    numericLabels = TRUE, 
    saveTOMs = T,
    verbose = 3
  )
  table(net$colors) 
  
}

if(T){
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "step3_genes_modules.Rdata")



rm(list = ls())
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
datTraits$group <- factor(datTraits$group, levels = c("CON", "ASD","BD","MDD","PTSD","SZ"))
datTraits$group
if(T){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group)
  
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes
  MEs <- orderMEs(MES0)
  
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("step4_Module-trait-relationship_heatmap2.pdf",
      width = 1.5*length(colnames(design)),
      height = 0.4*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = colorRampPalette(c("#105eb7", "white","#d7131a"))(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = "step4_design.Rdata")
}
if(T){
  y = as.data.frame(design[ ,"MDD"]);
  GS1=as.numeric(cor(y,datExpr,use="p"))
  GeneSignificance=abs(GS1)
  ModuleSignificance=tapply(GeneSignificance,moduleColors , mean , na.rm=T)
  pdf("step4-Module-single_trait-relationship_barplot_MDD.pdf",width = 8,height = 5)
  plotModuleSignificance(GeneSignificance,moduleColors)
  dev.off()
}

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
              legend = "") +stat_compare_means()
  }  
  colorNames <- names(MEs)
  pdf("step4_Module-trait-relationship_boxplot.pdf", width =16,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=4))
  dev.off()
}




levels(datTraits$group)
choose_group <- "SZ"

if(T){
  modNames <- substring(names(MEs), 3)
  
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  
  trait <- as.data.frame(design[,choose_group])
  geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  selectModule <- modNames
  pdf("step4_gene-Module-trait-significance_SZ.pdf",width=7, height=1.5*ncol(MEs))
  par(mfrow=c(ceiling(length(selectModule)/2),2))
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


rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

if(T){
  TOM=TOMsimilarityFromExpr(datExpr,power=power)
  dissTOM=1-TOM
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


if(T){
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  if(T){
    CON = as.data.frame(design[,1])
    names(CON) = "CON"
    MET = orderMEs(cbind(MEs, CON))
  }
  pdf("step5_module_cor_Eigengene-dendrogram.pdf",width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="",
                        marDendro = c(0,4,1,4),
                        marHeatmap = c(5,5,1,2),
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}

rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

OrgDb = "org.Hs.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"
table(moduleColors)



choose_module <- c("black","brown","cyan","green",
                   "grey","magenta","midnightblue","tan",
                   "turquoise","purple","red","salmon","yellow") 

if(T){
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  
  gene_module <- data.frame(gene=colnames(datExpr),
                            module=moduleColors)
  write.csv(gene_module,file = "step6_gene_moduleColors.csv",row.names = F, quote = F)
  tmp <- bitr(gene_module$gene,fromType = genetype,
              toType = "ENTREZID",
              OrgDb = OrgDb )
  gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  formula_res <- compareCluster(
    ENTREZID~module,
    data = choose_gene_module_entrz,
    fun = "enrichGO",
    OrgDb = OrgDb,
    ont = "BP",  #One of "BP", "MF", and "CC"  or "ALL"
    pAdjustMethod = "BH",   #指定多重假设检验矫正的方法,一般选择 "BH" 或 "fdr"，BH较严格，fdr较温和（计算的q小些）
    pvalueCutoff = 0.05
  )
  
  lineage1_ego <- simplify(
    formula_res,
    cutoff=0.5,
    by="p.adjust",
    select_fun=min
  )
  save(gene_module, formula_res, lineage1_ego, file="step6_module_GO_term.Rdata")
  write.csv(lineage1_ego@compareClusterResult,
            file="step6_module_GO_term.csv")
  dotp <- dotplot(lineage1_ego,
                  showCategory=10,
                  includeAll = TRUE,
                  label_format=90)
  ggsave(dotp,filename= "step6_module_GO_term.pdf", #device = cairo_pdf,
         width = 16,
         height = 15)
}

rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step3_genes_modules.Rdata")
table(moduleColors)


module = "brown"
if(T){
  dat=datExpr[,moduleColors==module]
  library(pheatmap)
  n=t(scale(dat))
  
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

rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")
module = "black" 
module = "magenta" 
module = "brown" 
module = "tan" 
module = "red" 

if(T){
  gene <- colnames(datExpr)
  inModule <- moduleColors==module
  modgene <- gene[inModule]
  write.table(modgene,paste0("step8_",module,"_modgene.csv"))
  
  TOM <- TOMsimilarityFromExpr(datExpr,power=power)
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modgene,modgene)
  
  nTop = 100#30
  IMConn = softConnectivity(datExpr[, modgene])
  top = (rank(-IMConn) <= nTop)
  filter_modTOM <- modTOM[top, top]
  
  vis <- exportNetworkToVisANT(filter_modTOM,
                               file = paste("step8_visANTinput-",module,".txt",sep = ""),
                               weighted = T,threshold = 0)
  cyt <- exportNetworkToCytoscape(filter_modTOM,
                                  edgeFile = paste("step8_CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("step8_CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                  weighted = TRUE,
                                  threshold = 0.15,
                                  nodeNames = modgene[top],
                                  nodeAttr = moduleColors[inModule][top])
}


modulegene <- read.csv("~/精神疾病/data/388/Seurat/data/WGCNA/模块：15个/step6_gene_moduleColors.csv")
modulegene <- subset(modulegene,subset = module %in% "brown") %>% pull(gene)
library(ggplot2)
library(ggpubr)
library(DOSE)
library(enrichplot)
library(clusterProfiler)
modulegene <- read.csv("~/精神疾病/data/388/Seurat/data/WGCNA/模块：15个/step6_gene_moduleColors.csv")
green <- subset(modulegene,subset = module %in% "brown") %>% pull(gene)

library(org.Hs.eg.db)
bp <-
  enrichGO(
    green,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term <- bp@result

library(ggplot2)
go <- c(
  "Golgi vesicle transport",#高尔基囊泡运输
  "RNA splicing",#RNA剪接
  "RNA localization",#RNA定位
  "histone modification",#组蛋白修饰
  "nucleic acid transport",#核酸转运
  "regulation of mRNA processing",#mRNA加工的调节
  "internal protein amino acid acetylation",#内部蛋白质氨基酸乙酰化
  "retrograde vesicle-mediated transport, Golgi to endoplasmic reticulum",
  "regulation of response to DNA damage stimulus",#DNA损伤刺激反应的调节
  "macroautophagy",#自噬
  "negative regulation of gene expression, epigenetic"#基因表达的负调控，表观遗传
  
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
pdf("~/精神疾病/data/388/Seurat/data/WGCNA/模块：15个/step9_GO_brown.pdf",height = 4,width = 6)
p <- ggplot(data = term,
            aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#CC6633",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=3.5,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
print(p)
dev.off()


