
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description='normalize count matrix using edgeR func. (lib size + TMM factor)')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-s', '--sample', type='character', default="all",
                    help='tsv/txt file: pass filter colnames/samples, default: all (keep all colnames/samples), optional: txt path of valid samples, will change mat to the same order')
parser$add_argument('-g', '--group', type='character', default="FALSE",
                    help='group info of colnames/samples, default: FALSE, (consider as one group), optional: txt path of group, in same order of valid samples')
parser$add_argument('-m', '--method', type='character', default="TMM",
                    help='normalize method: "TMM","TMMwsp","RLE","upperquartile","none", default: TMM')
parser$add_argument('-f', '--filterByExpr', type='character', default="FALSE",
                    help='remove low expr records/genes, default: FALSE, optional: TRUE')
parser$add_argument('-p', '--priCount', type='integer', default=2,
                    help='priCount: int var for cpm(), default: 2')
parser$add_argument('--geneLength', type='character', default="FALSE",
                    help='if provide geneLength txt path, then cal rpkm mat as output instead of cpm, default: FALSE, optional: txt path of geneLength, in same order of valid genes')
parser$add_argument('-e', '--externalLibSize', type='character', default="FALSE",
                    help='to build DEGlist replace colSums, default: FALSE, optional: txt path of libsize, in same order of valid samples')
parser$add_argument('-l', '--log', type='character', default="FALSE",
                    help='whether to log cpm result, default: FALSE, optional: TRUE')
parser$add_argument('-o', '--outfile', type='character', required=TRUE,
                    help='output file path')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}

ma <- args$matrix
sa <- (args$sample)
gr <- (args$group)
fi <- toupper(args$filterByExpr)
me <- args$method
pr <- args$priCount
gl <- (args$geneLength)
ex <- (args$externalLibSize)
lo <- args$log
ou <- args$outfile






message(paste0("read input mat from ",ma))
count.matrix <- read.csv(ma,sep = "\t",header = TRUE, check.names = F,row.names = 1)

if(sa=="all"){
  message(paste0("keep all samples from input mat"))
  s <- colnames(count.matrix)
} else {
  message(paste0("only keep samples from path: ",sa))
  s <- read.delim(sa,header = F,stringsAsFactors = F)$V1
}
count.matrix <- count.matrix[,s,drop=F]

if(gr=="FALSE"){
  message("not using group info (all as 1)")
  g <- rep(1,length(s))
} else {
  message(paste0("using group info from path: ",gr))
  g <- read.delim(gr,header = F,stringsAsFactors = F)$V1
}
tissue.types <- factor(g,level=unique(g)) 

if (ex=="FALSE"){
  message("use rowSums as lib.size in DGEList")
  y <- edgeR::DGEList(counts=count.matrix,group=tissue.types)
}else{
  message("use provided path as lib.size in DGEList")
  libSize <- read.table(ex,header = F,stringsAsFactors = F)
  tmp <- libSize$V1
  libSize <- libSize$V2
  names(libSize) <- tmp
  y <- edgeR::DGEList(counts=count.matrix,group=tissue.types,lib.size=libSize[s])
}
print(paste0("dim: ",dim(y$counts)))
print(paste0("grp: ",table(tissue.types)))

if(fi=="TRUE"){
  message("filter low expr by group info")
  keep <- edgeR::filterByExpr(y,group = tissue.types)
  if (ex=="FALSE"){
    message("libSize not kept !")
    y <- y[keep, , keep.lib.sizes=FALSE]
  } else{
    message("libSize is kept !")
    y <- y[keep, , keep.lib.sizes=TRUE]
  }
}else if(fi=="FALSE"){
  message("do not filter low expr")
}


message(paste0("cal scaling factor using ",me))
y <- edgeR::calcNormFactors(y, method=me)
if(gl=='FALSE'){
  message(paste0("gene length not provided, will output cpm mat"))
  y.cpm <- as.data.frame(edgeR::cpm(y,normalized.lib.sizes = TRUE,
                                    log = lo, 
                                    prior.count = pr) 
                         )
} else {
  message(paste0("gene length provided, will output rpkm mat"))
  geneLen <- read.delim(gl,header = F,stringsAsFactors = F)$V1
  y.cpm <- as.data.frame(edgeR::rpkm(y,normalized.lib.sizes = TRUE,
                                     log = lo, 
                                     prior.count = pr, 
                                     gene.length = geneLen)
                         )
}




y.cpm <- cbind(rownames(y.cpm),y.cpm)
colnames(y.cpm)[1] <- "gene_id"

message(paste0("write file to ",ou))
write.table(y.cpm, file = ou, sep = "\t",col.names = T,row.names = F,quote = F)
