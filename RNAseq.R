##############################################################################
######################## loading the environment #############################
##############################################################################
BiocManager::install("edgeR")
#loading in files
counts <- as.matrix(read.table("RNA_seq/software/RNA_seq_analysis-master/htseq.count.txt",sep="\t", 
                               header =T ,row.names = 1))
group <- as.character(read.table("RNA_seq/software/RNA_seq_analysis-master/sample_list.txt")[,1])

#loading data into DGElist class
library(edgeR)
cds <- DGEList( counts, group = group )

##############################################################################
######################## transforming the data ###############################
##############################################################################
#removing genes that have low expression values
cds <- cds[rowSums(1e+06 * (cds$counts/expandAsMatrix(cds$samples
                                                      $lib.size, dim(cds)))> 1) >= 2, ]
cds <- calcNormFactors( cds )
#exporting high read files
CPM <- cpm(cds)
write.table(CPM,"RNA_seq/analysis/CPM.txt",sep="\t",quote=F)

##############################################################################
######################## stats and significance ##############################
##############################################################################
#calculating fold change
#################### no comparison to be had
cds <- estimateCommonDisp(cds)
de.poi <- exactTest( cds, pair = c( "control", "GD"))
resultsByFC.poi <- topTags( de.poi, n = nrow( de.poi$table ), sort.by
                            = "logFC" )$table

#identify and create differential gene file
de.gene.logFC <- resultsByFC.poi[abs(resultsByFC.poi[,"logFC"])>1 &
                                   resultsByFC.poi[,"FDR"]<0.01, ]
write.table(de.gene.logFC,"DE.gene.logFC.xls",sep="\t",quote=F)

#visiual of differential expression
png("figures/meet16_RNAseq/Volcano.png",5,5,units = "in",res=300)
plot(resultsByFC.poi[,1],-log10(resultsByFC.poi[,"FDR"]),xlab="logFC",ylab="-log(P-value)",pch=20,cex=0.3)
points(de.gene.logFC[,1],-log10(de.gene.logFC[,"FDR"]),col="red",pch=1,cex=0.3)
dev.off()

########
#heatmap of gene expredssion
library(gplots)
de.genes = rownames(de.gene.logFC)
rpkm = rpkm(cds, counts)
de.gene.RPKM.mtx = log10(rpkm[de.genes,]+0.1)
colfunc <- colorRampPalette(brewer.pal(9,"Blues"))
de.gene.RPKM.mtx[is.na(de.gene.RPKM.mtx)] = 0
png("figures/meet16_RNAseq/DE.gene.heatmap.png",5,7,units = "in",res=300) 
heatmap.2(de.gene.RPKM.mtx,col=colfunc(100),trace="none",labRow="",cexRow=1,cexCol=1.5,lhei=c(1,5), lwid=c(1, 3), density.info="none",margins =c(5,1),srtCol=45,key.title = "",key.xlab = "logRPKM")
dev.off()




