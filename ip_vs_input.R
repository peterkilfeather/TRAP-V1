library("DESeq2")
library("tximport")
library("GenomicFeatures")
library("biomaRt")
library("reshape")
library(reshape2)
library(ggplot2)
library(pheatmap)

wd <- "/Users/peterkilfeather/Desktop/pk_trap"
setwd(wd)

sample_metadata <- read.csv("sample_metadata_020619.csv", header = T)
sample_metadata_striatum_all <- sample_metadata[sample_metadata$region != 1,]
sample_metadata_striatum_all$collection_day[is.na(sample_metadata_striatum_all$collection_day)] <- 10
sample_metadata_striatum_all$collection_day <- factor(sample_metadata_striatum_all$collection_day)
sample_metadata_striatum_all$age_months <- factor(sample_metadata_striatum_all$age_months)
sample_metadata_striatum_all$ip <- factor(sample_metadata_striatum_all$ip)
sample_metadata_striatum_all$region <- factor(sample_metadata_striatum_all$region)

files_striatum_all = file.path(wd, "quant_pc", sample_metadata$sample_code, "abundance.h5")[sample_metadata$region != 1]
names(files_striatum_all) <- sample_metadata$sample_name[sample_metadata$region != 1]

txdb <- makeTxDbFromGFF("Mus_musculus.GRCm38.96.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

txi_ip_vs_input <- tximport(files_striatum_all, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)

dds_ip_vs_input <- DESeqDataSetFromTximport(txi_ip_vs_input, colData = sample_metadata_striatum_all, design = ~ age_months + ip)

dim(dds_ip_vs_input)


hist(as.numeric(rowSums(counts(dds_ip_vs_input))),
     breaks = 100, main = "Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(counts(dds_ip_vs_input)))),col="red")

# ribosomal <- rowMeans(counts(dds_ip_vs_input)) > 150000
# sum(ribosomal, na.rm = T)
# high_count <- dds_ip_vs_input[ribosomal,]
# dim(dds_IPvsT)

keep_feature <- rowSums(counts(dds_ip_vs_input)) > 0
dds_ip_vs_input <- dds_ip_vs_input[keep_feature, ]
dim(dds_ip_vs_input)

# density plot of raw read counts (log10)
counts_ip_vs_input <- counts(dds_ip_vs_input)
logcounts <- log(counts_ip_vs_input[,1]+1,10) 
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.45),xlab="Raw read counts per gene (log10)", ylab="Density")
for (s in 2:ncol(counts_ip_vs_input)){
  logcounts <- log(counts_ip_vs_input[,s],10) 
  d <- density(logcounts)
  lines(d)
}

counts_ip_vs_input_log = log(counts_ip_vs_input+1)
hist(as.numeric(rowSums(counts_ip_vs_input_log)),
     breaks = 100, main = "log Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(counts_ip_vs_input_log))),col="red")

select = order(rowMeans(counts_ip_vs_input), decreasing=TRUE)[1:100]
highexprgenes_counts <- counts_ip_vs_input[select,]
highexprgenes_logcounts <- counts_ip_vs_input_log[select,]

hist(as.numeric(rowSums(highexprgenes_logcounts)),
     breaks = 100, main = "log Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(highexprgenes_logcounts))),col="red")

par(mfrow=c(2,2)) ## if you want to have multiple plots on the same window.

# heatmap with sample name on X-axis
heatmap(as.matrix(highexprgenes_logcounts), col=topo.colors(50), margin=c(10,6))
# heatmap with condition group as labels
colnames(highexprgenes_logcounts)<- sample_metadata_striatum_all$sample_name
# plot
heatmap(as.matrix(highexprgenes_logcounts), col = topo.colors(50), margin=c(10,6))

data_for_PCA <- t(highexprgenes_logcounts)
dim(data_for_PCA)

############### Back to slides #################
data_for_PCA_scaled = scale(data_for_PCA,center = TRUE, scale = FALSE)

# using prcomp
bpc = prcomp(data_for_PCA, center=TRUE, scale=FALSE) ## stats package 
beigenvalues = bpc$sdev^2 # eigenvalues
beigenvectors = bpc$rotation # eigenvectors
par(mfrow=c(2,1))
plot(beigenvalues/sum(beigenvalues) * 100,xlab="PC",ylab="% Variance explained") 
plot(cumsum(beigenvalues)/sum(beigenvalues) * 100, xlab="PC",ylab="Cumulative % Variance explained")

plot(bpc$x[,1]/(bpc$sdev[1]*sqrt(12)),bpc$x[,2]/(bpc$sdev[2]*sqrt(12)), xlab="PC1",ylab="PC2", main="PC prcomp")
library(ggfortify)
autoplot(bpc, main = "PC prcomp, autoplot")

library(pcaExplorer)
# extract genes with the highest/lowest loadings per PC
hi_loadings(bpc, whichpc = 1, topN = 10, exprTable = NULL,
            annotation = NULL, title = "Top/bottom loadings - ")

hi_loadings(bpc, whichpc = 2, topN = 10, exprTable = NULL,
            annotation = NULL, title = "Top/bottom loadings - ")

s = data_for_PCA[,which(colnames(data_for_PCA)=="ENSMUSG00000029580")]
s2 = cbind(s,sample_metadata_striatum_all)
ggplot(s2, aes(ip,s, fill=as.factor(region)))+ geom_boxplot(position="dodge") +geom_point(alpha=0.6, aes(group=age_months), data=s2, position = position_dodge(width=0.75))+ylab("log Expression")+xlab("")+ggtitle("Gene 213742")

# doing the above with custom code
whichpc = 1
topN = 10
geneloadings_sorted <- sort(bpc$rotation[, whichpc])
geneloadings_extreme <- c(tail(geneloadings_sorted, topN), head(geneloadings_sorted, topN))
barplot(geneloadings_extreme, las = 2, col = c(rep("steelblue", topN), rep("coral", topN)), main = paste0("PC", whichpc," loadings"))

# plotting a heatmap of expression values for the top topN genes
s = data_for_PCA[,which(colnames(data_for_PCA)%in%names(geneloadings_extreme))]
s2 = cbind(s,sample_metadata_striatum_all)
data_narrow <- melt(s2, id =names(sample_metadata_striatum_all))
ggplot(data_narrow, aes(ip,variable))+geom_tile(aes(fill=value))+facet_grid(.~region)+ylab("Gene")+scale_fill_gradientn(colours = rainbow(20))  #+ scale_fill_gradient2( low = "blue", high = "red", na.value="black", name = "" )
pheatmap(s)

# doing the same for PC2
whichpc = 2
topN = 10
geneloadings_sorted <- sort(bpc$rotation[, whichpc])
geneloadings_extreme <- c(tail(geneloadings_sorted, topN), head(geneloadings_sorted, topN))
barplot(geneloadings_extreme, las = 2, col = c(rep("steelblue", topN), rep("coral", topN)), main = paste0("PC", whichpc," loadings"))
s = data_for_PCA[,which(colnames(data_for_PCA)%in%names(geneloadings_extreme))]
s2 = cbind(s,sample_metadata_striatum_all)
data_narrow <- melt(s2, id =names(sample_metadata_striatum_all))
ggplot(data_narrow, aes(ip,variable))+geom_tile(aes(fill=value))+facet_grid(.~region)+ylab("Gene")+scale_fill_gradientn(colours = rainbow(20))  #+ scale_fill_gradient2( low = "blue", high = "red", na.value="black", name = "" )
pheatmap(s)

###
# plot gene expression on top of PCA plot
s = data_for_PCA[,which(colnames(data_for_PCA)=="ENSMUSG00000044308")]
logExpression = as.numeric(s)
shape = paste(sample_metadata_striatum_all$region,sep="-")
ggplot(bpc$x, aes(PC1,PC2,colour=logExpression)) +geom_point(aes(shape=shape)) +ggtitle("Ubr3 expression")

# other ways to perform tsne - slow so not run here for the scRNA-seq data
library("Rtsne")
matrix <- as.matrix(data_for_PCA)
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(matrix, dims = 2, perplexity = 9) # Run TSNE
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=sample_metadata_striatum_all$region, pch=as.numeric(sample_metadata_striatum_all$ip))

dat = tsne_out$Y
colnames(dat) = c("tSNE1","tSNE2")
ggplot(dat, aes(tSNE1,tSNE2))+geom_point() +ggtitle("tSNE")

s = data_for_PCA[,which(colnames(data_for_PCA)=="ENSMUSG00000044308")]
logExpression = as.numeric(s)
shape = paste(sample_metadata_striatum_all$ip,sample_metadata_striatum_all$region,sep="-")
ggplot(dat, aes(tSNE1,tSNE2,colour=logExpression)) +geom_point(aes(shape=shape)) +ggtitle("Ubr3 expression")

s = data_for_PCA[,which(colnames(data_for_PCA)=="ENSMUSG00000044308")]
logExpression = as.numeric(s)
s2 = cbind(s,sample_metadata_striatum_all)
s2$cluster = paste(s2$ip,":",s2$region, sep="")
ggplot(s2,aes(cluster, s,fill=cluster))+geom_boxplot()+geom_point()

# kmeans clustering https://www.r-bloggers.com/k-means-clustering-in-r/
Cluster <- kmeans(data_for_PCA_scaled[,1:2], 3, nstart = 20)
Cluster$cluster <- as.factor(Cluster$cluster)
# get cluster means 
aggregate(data_for_PCA[,1:2],by=list(Cluster$cluster),FUN=mean)

ggplot(data_for_PCA[,1:2], aes(`ENSMUSG00000029580`,`ENSMUSG00000015656`, color = Cluster$cluster, shape=shape)) + geom_point()

par(mfrow=c(1,1))
measure= NULL
for (k in 1:11){
  Clust <- kmeans(data_for_PCA[,1:2],k,nstart=12)
  measure = c(measure,Clust$betweenss/Clust$totss)
}
plot(measure, xlab="k", ylab="between_SS / total_SS", type="b", pch=19)

#clusters 
Cluster$cluster
# their centers
Cluster$centers
# numbers of points per cluster
Cluster$size

library(mclust)
d_clust <- Mclust(as.matrix(data_for_PCA_scaled), G=1:9, 
                  modelNames = mclust.options("emModelNames"))
d_clust$BIC
plot(d_clust)

library(NbClust)
nb <- NbClust(data_for_PCA_scaled, diss=NULL, distance = "euclidean", 
              min.nc=2, max.nc=5, method = "kmeans", 
              index = "all", alphaBeale = 0.1)
hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))


quantile(rowMeans(counts(dds_ip_vs_input)), seq(0, 1, 0.05))
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_ip_vs_input)[["cooks"]]), range=0, las=2)

keep <- rowMeans(counts(dds_ip_vs_input)) >= 100 # must have 100 counts  per sample for consideration
# keep <- rowSums((counts(dds_IPvsT)) >= 10 ) > 10 
sum(keep, na.rm = T)
# sum(keep1, na.rm = T)
# sum(keep2, na.rm = T)
dds_ip_vs_input <- dds_ip_vs_input[keep,]
dim(dds_ip_vs_input)

dds_ip_vs_input <- DESeq(dds_ip_vs_input, minReplicatesForReplace = Inf)

plotDispEsts(dds_ip_vs_input)
abline(h=0.1,col="green")

res_ip_vs_input <- results(dds_ip_vs_input, alpha = 0.001)
summary(res_ip_vs_input)

res_ip_vs_input$ensembl_gene_id = rownames(res_ip_vs_input)
res_ip_vs_input_enriched <- subset(as.data.frame(res_ip_vs_input), padj < 0.001)
res_ip_vs_input_enriched <- res_ip_vs_input[order(-res_ip_vs_input$baseMean),]
res_ip_vs_input_enriched <- merge(as.data.frame(res_ip_vs_input_enriched), anno, by=c("ensembl_gene_id"))
write.csv(res_ip_vs_input_enriched, file = "striatum_enriched.csv")


mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
anno <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), filters = "ensembl_gene_id", values = res_ip_vs_input$ensembl, mart = mart)

fasn_coverage <- read.table("fasn.coverage", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

fasn_coverage <- rename(fasn_coverage,c(V1="Chr", V2="locus", V3="depth")) # renames the header
plot(fasn_coverage$locus, fasn_coverage$depth)
library(lattice, pos=10) 
xyplot(depth ~ locus, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), data=fasn_coverage, main="sequence coverage across Fasn locus")


