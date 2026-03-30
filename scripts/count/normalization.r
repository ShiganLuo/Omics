library(edgeR)

TMM = function(counts) {
    count_data = read.csv(counts,header=TRUE,sep="\t",row.names=1)
    count_matrix <- as.matrix(count_data)
    
}