library("DESeq2")
library("tximport")
library("GenomicFeatures")
library("biomaRt")
library("reshape")

wd <- "/Users/peterkilfeather/Desktop/pk_trap"
setwd(wd)

sample_metadata <- read.csv("sample_metadata_020619.csv", header = T)
sample_metadata_striatum_all <- sample_metadata[sample_metadata$region != 1,]
sample_metadata_striatum_all$collection_day[is.na(sample_metadata_striatum_all$collection_day)] <- 10
sample_metadata_striatum_all$collection_day <- factor(sample_metadata_striatum_all$collection_day)
sample_metadata_striatum_all$age_months <- factor(sample_metadata_striatum_all$age_months)
sample_metadata_striatum_all$ip <- factor(sample_metadata_striatum_all$ip)

files_striatum_all = file.path(wd, "quant_pc", sample_metadata$sample_code, "abundance.h5")[sample_metadata$region != 1]
names(files_striatum_all) <- sample_metadata$sample_name[sample_metadata$region != 1]

txdb <- makeTxDbFromGFF("Mus_musculus.GRCm38.96.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

txi_ip_vs_input <- tximport(files_striatum_all, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)

dds_ip_vs_input <- DESeqDataSetFromTximport(txi_ip_vs_input, colData = sample_metadata_striatum_all, design = ~ age_months + ip)

dim(dds_ip_vs_input)

ribosomal <- rowMeans(counts(dds_ip_vs_input)) > 150000
sum(ribosomal, na.rm = T)
high_count <- dds_ip_vs_input[ribosomal,]
dim(dds_IPvsT)

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


