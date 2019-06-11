library("DESeq2")
library("tximport")
library("GenomicFeatures")
library("biomaRt")
library("reshape")
library(reshape2)
library(ggplot2)
library(pheatmap)
library(scater)
library("RColorBrewer")
library(dplyr)
library("org.Mm.eg.db")

require(scatterplot3d)

mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
anno <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), filters = "ensembl_gene_id", values = rownames(counts), mart = mart)

sample_metadata <- read.csv("sample_metadata_020619.csv", header = T)
sample_metadata_striatum_all <- sample_metadata[sample_metadata$region != 1,] 
sample_metadata_striatum_all$collection_day[is.na(sample_metadata_striatum_all$collection_day)] <- 10
sample_metadata_striatum_all$collection_day <- factor(sample_metadata_striatum_all$collection_day)
sample_metadata_striatum_all$age_months <- factor(sample_metadata_striatum_all$age_months)
sample_metadata_striatum_all$ip <- factor(sample_metadata_striatum_all$ip)
sample_metadata_striatum_all$region <- factor(sample_metadata_striatum_all$region)

files_striatum_all = file.path("/zfs/analysis/PK_TRAP/", "quant_pc", sample_metadata$sample_code, "abundance.h5")[sample_metadata$region != 1]
names(files_striatum_all) <- sample_metadata$sample_name[sample_metadata$region != 1]

txdb <- makeTxDbFromGFF("Mus_musculus.GRCm38.96.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

txi_ip_vs_input <- tximport(files_striatum_all, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)

dds_ip_vs_input <- DESeqDataSetFromTximport(txi_ip_vs_input, colData = sample_metadata_striatum_all, design = ~ age_months + region + ip)
dds_ip_vs_input <- DESeq(dds_ip_vs_input, minReplicatesForReplace = Inf)
dim(dds_ip_vs_input)

#filter for min 10 counts in min 2 samples
keep_feature <- rowSums(counts(dds_ip_vs_input) >= 10) >=2
dds_ip_vs_input <- dds_ip_vs_input[keep_feature, ]
dim(dds_ip_vs_input)

#variance stabilised counts: Note that blind is set to False
vsd_ip_vs_input <- vst(dds_ip_vs_input, blind=F)

#Choose
counts <- counts(dds_ip_vs_input, normalized = T)
counts <- assay(vsd_ip_vs_input)

# density plot of raw read counts (log10)
logcounts <- log(counts[,1]+1,10) 
d <- density(logcounts)
plot(d,xlim=c(0,8),main="",ylim=c(0,.45),xlab="Raw read counts per gene (log10)", ylab="Density")
for (s in 2:ncol(counts)){
  logcounts <- log(counts[,s]+1,10) 
  d <- density(logcounts)
  lines(d)
}

# density plot of raw read counts (vsd)
logcounts <- log10(counts)
d <- density(logcounts)
plot(d,xlim=c(0,2),main="",ylim=c(0,30),xlab="Raw read counts per gene (log10)", ylab="Density")
for (s in 2:ncol(counts)){
  logcounts <- log(counts[,s],10) 
  d <- density(logcounts)
  lines(d)
}

#expression sum per gene
logcounts = log10(counts+1)
hist(as.numeric(rowSums(logcounts)),
     breaks = 100, main = "log Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(counts_log))),col="red")

#sample distances
sampleDists_ip_vs_input <- dist(t(assay(vsd_ip_vs_input)))
sampleDistMatrix_ip_vs_input <- as.matrix(sampleDists_ip_vs_input)
rownames(sampleDistMatrix_ip_vs_input) <- paste(vsd_ip_vs_input$sample_name)
colnames(sampleDistMatrix_ip_vs_input) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_ip_vs_input,
         clustering_distance_rows=sampleDists_ip_vs_input,
         clustering_distance_cols=sampleDists_ip_vs_input,
         col=colors)

##PCA with vsd for ip vs input
#calculate the variance for each gene
feature_variance <- rowVars(assay(vsd_ip_vs_input))

## select the ntop genes by variance
select <- order(feature_variance, decreasing=TRUE)[seq_len(min(500, length(feature_variance)))]
#select = order(rowMeans(counts), decreasing=TRUE)[1:100]
# highexprgenes_counts <- counts[select,]
# highexprgenes_logcounts <- counts_log[select,]

highvargenes_counts <- counts[select,]
highvargenes_logcounts <- logcounts[select,]

hist(as.numeric(rowSums(highvargenes_logcounts)),
     breaks = 100, main = "log Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(highvargenes_logcounts))),col="red")

par(mfrow=c(2,2)) ## if you want to have multiple plots on the same window.

# heatmap with sample name on X-axis
heatmap(as.matrix(highvargenes_logcounts), col=topo.colors(50), margin=c(10,6))
# heatmap with condition group as labels
colnames(highvargenes_logcounts)<- sample_metadata_striatum_all$sample_name
# plot
heatmap(as.matrix(highvargenes_logcounts), col = topo.colors(50), margin=c(10,6))

###OBDS PCA
data_for_PCA <- t(highvargenes_logcounts)
dim(data_for_PCA)

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

s = data_for_PCA[,which(colnames(data_for_PCA)=="ENSMUSG00000067879")]
s2 = cbind(s,sample_metadata_striatum_all)
ggplot(s2, aes(ip,s, fill=as.factor(region)))+ geom_boxplot(position="dodge") + geom_point(alpha=0.6, aes(group=region), data=s2, position = position_dodge(width=0.75))+ylab("log Expression")+xlab("")+ggtitle("Vexin")

###Return to finding IP enriched genes

plotDispEsts(dds_ip_vs_input)
abline(h=0.1,col="green")

res_ip_vs_input <- results(dds_ip_vs_input, alpha = 0.001)
summary(res_ip_vs_input)

res_ip_vs_input$ensembl_gene_id = rownames(res_ip_vs_input)
res_ip_vs_input_enriched <- subset(as.data.frame(res_ip_vs_input), padj < 0.001 & log2FoldChange > 0.25)
res_ip_vs_input_enriched <- res_ip_vs_input_enriched[order(res_ip_vs_input_enriched$padj),]
res_ip_vs_input_enriched <- merge(as.data.frame(res_ip_vs_input_enriched), anno, by=c("ensembl_gene_id"))
write.csv(res_ip_vs_input_enriched, file = "striatum_enriched.csv")

plotCounts(dds_ip_vs_input, gene="ENSMUSG00000019943", intgroup="region")

striatal_genes <- rownames(res_ip_vs_input_enriched)

#just IP samples
txi_ip <- tximport(files_striatum_all[1:20], type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
sample_metadata <- sample_metadata_striatum_all[1:20,]

dds_ip <- DESeqDataSetFromTximport(txi_ip, colData = sample_metadata, design = ~ collection_day + age_months + region)
#subset for striatal enriched genes
striatal_keep <- rownames(dds_ip) %in% striatal_genes
dds_ip <- dds_ip[striatal_keep, ]
dim(dds_ip)

#filter min 10 counts in min 2 samples
keep_feature <- rowSums(counts(dds_ip) >= 10) >=2
dds_ip <- dds_ip[keep_feature, ]
dim(dds_ip)

dds_ip <- DESeq(dds_ip, minReplicatesForReplace = Inf)

res_ip_region <- results(dds_ip)
summary(res_ip_region)

# df <- data.frame(Sample=factor(c(1:10)), Group=c(rep("A",5), rep("B",5)))
# mm <- model.matrix(~Sample+Group, df)
# mm

res_ip_region$symbol <- mapIds(x=org.Mm.eg.db, keys = rownames(res_ip_region), keytype="ENSEMBL", column="SYMBOL")
res_ip_region$ensembl <- rownames(res_ip_region)
res_ip_region <- as.data.frame(res_ip_region) %>% filter(symbol != 'NA') %>% arrange(padj, -log2FoldChange)

res_ip_ventral_enriched <- subset(as.data.frame(res_ip_region), padj < 0.001 & log2FoldChange > 0.25)
res_ip_dorsal_enriched <- subset(as.data.frame(res_ip_region), padj < 0.001 & log2FoldChange < 0.25)
  
head(res_ip_dorsal_enriched, 20)
head(res_ip_ventral_enriched, 20)

res_ip_dorsal_enriched <- merge(as.data.frame(res_ip_dorsal_enriched), anno, by=c("ensembl_gene_id"))

res_ip_dorsal_enriched <- left_join(res_ip_dorsal_enriched, anno, by = c("ensembl" = "ensembl_gene_id"))
res_ip_ventral_enriched <- left_join(res_ip_ventral_enriched, anno, by = c("ensembl" = "ensembl_gene_id"))

write.csv(res_ip_dorsal_enriched, file = "dorsal_striatum_enriched.csv")
write.csv(res_ip_ventral_enriched, file = "ventral_striatum_enriched.csv")


###Generic PCA

## perform a PCA on the data in assay(x) for the selected genes
pca_striatum <- prcomp(t(assay(vsd_ip_vs_input)[select,]))

## the contribution to the total variance for each component
percentVar_striatum <- pca_striatum$sdev^2 / sum( pca_striatum$sdev^2 )

##plot the "percentVar"
scree_plot=data.frame(percentVar_striatum)
scree_plot[,2]<- c(1:28)

colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")



###PCA Kevin Blighe
pca_striatum <- prcomp(t(assay(vsd_ip_vs_input)))
summary(pca_striatum)
pca_striatum_proportionvariances <- ((pca_striatum$sdev^2) / (sum(pca_striatum$sdev^2)))*100
barplot(pca_striatum_proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca_striatum$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(pca_striatum$x[,1:5], col="blue", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
pairs(pca_striatum$x[,6:10], col="blue", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)

par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)

#Plots scatter plot for PC 1 and 2
plot(pca_striatum$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(pca_striatum_proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(pca_striatum_proportionvariances[2], 2), "%"))
points(pca_striatum$x, col="black", pch=16, cex=1)

#Plots scatter plot for PC 1 and 3
plot(pca_striatum$x[,1], pca_striatum$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(pca_striatum_proportionvariances[1], 2), "%"), ylab=paste("PC3, ", round(pca_striatum_proportionvariances[3], 2), "%"))
points(pca_striatum$x[,1], pca_striatum$x[,3], col="black", pch=16, cex=1)

#Plots scatter plot for PC 2 and 3
plot(pca_striatum$x[,2], pca_striatum$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC2, ", round(pca_striatum_proportionvariances[2], 2), "%"), ylab=paste("PC3, ", round(pca_striatum_proportionvariances[3], 2), "%"))
points(pca_striatum$x[,2], pca_striatum$x[,3], col="black", pch=16, cex=1)

par(mar=c(4,4,4,4), cex=1.0, cex.main=0.8, cex.axis=0.8)

scatterplot3d(pca_striatum$x[,1:3], angle=-40, main="", color="black", pch=17, xlab=paste("PC1, ", round(pca_striatum_proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(pca_striatum_proportionvariances[2], 2), "%"), zlab=paste("PC3, ", round(pca_striatum_proportionvariances[3], 2), "%"), grid=FALSE, box=FALSE)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca_striatum$x[,1:3], grid = c("xy", "xz", "yz"))

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca_striatum$x[,1:3], grid = c("xy", "xz", "yz"))

###End of PCA Kevin Blighe

hist(as.numeric(rowSums(counts(dds_ip_vs_input))),
     breaks = 100, main = "Expression sum per gene",
     xlab = "Sum expression")
abline(v=median(as.numeric(rowSums(counts(dds_ip_vs_input)))),col="red")



#scater
sce_all <- SingleCellExperiment(assays = list(counts = counts(dds_ip_vs_input)), colData = sample_metadata_striatum_all)
sce_all

plot(quantile(rowMeans(counts(dds_ip_vs_input)), seq(0, 1, 0.05))[1:20])
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_ip_vs_input)[["cooks"]]), range=0, las=2)

keep <- rowMeans(counts(dds_ip_vs_input)) >= 100 # must have 100 counts  per sample for consideration
# keep <- rowSums((counts(dds_IPvsT)) >= 10 ) > 10 
sum(keep, na.rm = T)
# sum(keep1, na.rm = T)
# sum(keep2, na.rm = T)
dds_ip_vs_input <- dds_ip_vs_input[keep,]
dim(dds_ip_vs_input)




fasn_coverage <- read.table("fasn.coverage", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

fasn_coverage <- rename(fasn_coverage,c(V1="Chr", V2="locus", V3="depth")) # renames the header
plot(fasn_coverage$locus, fasn_coverage$depth)
library(lattice, pos=10) 
xyplot(depth ~ locus, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), data=fasn_coverage, main="sequence coverage across Fasn locus")







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

# Library does not load, currently!
# library(NbClust)
# nb <- NbClust(data_for_PCA_scaled, diss=NULL, distance = "euclidean", 
#               min.nc=2, max.nc=5, method = "kmeans", 
#               index = "all", alphaBeale = 0.1)
# hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

#hierarchical clustering

d_pca <- dist(data_for_PCA) # euclidean by default
hc_pca <- hclust(d_pca, method = "complete") ## hierarchical cluster analysis
plot(hc_pca)
groups <- cutree(hc_pca,k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(hc_pca, k=5, border="red")

library(pheatmap)  
pheatmap(data_for_PCA[,1:10])

pheatmap(data_for_PCA)

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(t(data_for_PCA), method.hclust="ward",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

# Model Based Clustering
library(mclust)
fit <- Mclust(t(data_for_PCA))
plot(fit) # plot results 
summary(fit) # display the best model







