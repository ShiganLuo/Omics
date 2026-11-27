library(DiffBind)

#1. Load samples into DBA object
##Note that FACTOR column has cell line + treatment for per cell analysis of DARs
dba.sample=DiffBind::dba(sampleSheet="samplesheet_GSC_ATAC.v4.csv")
dba.sample=DiffBind::dba.count( dba.sample,
                                score=DBA_SCORE_READS,
                                fragmentSize = 1,
                                summits=250,
                                mapQCth = 10,
                                bParallel = T
                                )
gr.raw_peak=dba.peakset(dba.sample, bRetrieve = T ); 
df.raw_peak=as(gr.raw_peak, "data.frame")

df.ATAC.HouseKG=read.table("ATAC.housekeeping.bed")
colnames(df.ATAC.HouseKG)=c("seqnames", "start", "end", "peak_number", "score", "strand")
mat=df.raw_peak[ rownames(df.raw_peak) %in% df.ATAC.HouseKG$peak_number, 6:33]
vecGeoMean=apply( mat, 1, FUN=function(a_x){ return( exp( mean(log( a_x+1) ) ) ); })
mat.norm=sweep(mat, 1, vecGeoMean, FUN="/")
vecSizeFactor=apply( mat.norm, 2, FUN=function(a_x){ return( median(a_x) ); })
write.table(as.data.frame(vecSizeFactor), "sizefactor.txt", quote=F, row.names=T, col.names=T, sep="\t")



