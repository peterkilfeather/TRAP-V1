enrichment_cutoff = 3

# assay(vsd_ds_ip_vs_mb_ip)

#Select genes ip/input
axon_ip <- res_output_axon_ip_vs_axon_ub %>%
  filter(padj < 0.01)
ds_ip <- res_output_ds_ip_vs_ds_ub %>%
  filter(padj < 0.01)
vs_ip <- res_output_vs_ip_vs_vs_ub %>%
  filter(padj < 0.01)
mb_ip <- res_output_mb_ip_vs_mb_ub %>%
  filter(padj < 0.01)

#Select genes ds/vs/mb from ip/input

ds_deseq <- res_output_ds_ip_vs_mb_ip %>%
  filter(gene %in% ds_ip$gene & log2FoldChange > 3 & padj < 0.01)
vs_deseq <- res_output_vs_ip_vs_mb_ip %>%
  filter(gene %in% vs_ip$gene & log2FoldChange > 3 & padj < 0.01)
axon_deseq <- ds_deseq %>%
  filter(gene %in% vs_deseq$gene)
mb_deseq0 <- res_output_ds_ip_vs_mb_ip %>%
  filter(gene %in% mb_ip$gene & log2FoldChange < -3 & padj < 0.01)
mb_deseq1 <- res_output_vs_ip_vs_mb_ip %>%
  filter(gene %in% mb_ip$gene & log2FoldChange < -3 & padj < 0.01)
mb_deseq <- mb_deseq0 %>% filter(gene %in% mb_deseq1$gene)
                                 

vsd_counts <- assay(vsd_axon_ip_vs_mb_ip)
log_counts <- log2(counts_axon_ip_vs_mb_ip + 1)

sample_vs_sample <- log_counts %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>% 
  as_tibble %>% 
  # filter(gene %in% axon_ip$gene) %>%
  gather(sample, counts, -gene) %>%
  inner_join(sample_metadata, by = c("sample" = "sample_name")) %>%
  dplyr::select(gene, counts, region) %>%
  group_by(gene, region) %>%
  summarise(counts = mean(counts)) %>%
  pivot_wider(names_from = region, values_from = counts) %>%
  pivot_longer(cols = c('ds', 'vs'), names_to = 'striatum', values_to = 'striatum_counts') %>%
  select(mb_counts = mb)
  
  # rowwise() %>%
  # mutate(axon = mean(c(ds, vs), na.rm = T)) %>%
  # mutate(ds_vs_diff = ifelse(ds - vs >= 1, 'large', 'small')) %>%
  # mutate(ds_ratio = ds - mb) %>%
  # mutate(vs_ratio = vs - mb) %>%
  # mutate(specific = ifelse(ds_ratio < enrichment_cutoff & vs_ratio < enrichment_cutoff, 'none', ifelse(ds_ratio > enrichment_cutoff & vs_ratio > enrichment_cutoff, 'both', ifelse(ds_ratio > enrichment_cutoff, 'ds', 'vs')))) %>%
  # mutate(axon_specific = ifelse(axon - mb >= 3, 'axon', ifelse(mb - axon >= 3, 'mb', 'neither'))) %>%
  # mutate(specific_deseq = ifelse(gene %in% axon_deseq$gene, 'axon', ifelse(gene %in% mb_deseq$gene, 'mb', 'neither')))

p <- ggplot(sample_vs_sample, aes(mb, axon))
p + geom_point(alpha = 1/2, aes(colour = factor(ds_vs_diff))) +
  facet_wrap(~ specific_deseq) +
  geom_abline(intercept = 3, slope = 1, linetype="dotted") +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = -3, slope = 1, linetype="dotted") +
  theme_bw() +
  scale_color_npg() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab(expression('Log'[2]*' Midbrain IP Counts')) +
  ylab(expression('Log'[2]*' Axon IP Counts')) +
  labs(colour = 'Selective expression')


p <- ggplot(res_output_axon_ip_vs_mb_ip, aes(gene, lfcSE))
p + geom_point(alpha = 1/2, aes(colour = factor(padj < 0.01))) +
  theme_bw() +
  scale_color_npg()  +
  coord_cartesian(ylim = c(0, 6)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_abline(intercept = 0.58, slope = 0, linetype="dotted")

+
  geom_abline(intercept = 3, slope = 1, linetype="dotted") +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = -3, slope = 1, linetype="dotted") +
  theme_bw() +
  scale_color_npg() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab(expression('Log'[2]*' Midbrain IP Counts')) +
  ylab(expression('Log'[2]*' Axon IP Counts')) +
  labs(colour = 'Selective expression')

plot(log2( 1 + counts(dds_axon_ip_vs_mb_ip, normalized=TRUE)[ , c(7, 4)] ),
     pch=16, cex=0.3)



axon_enriched <- sample_vs_sample %>%
  filter(specific %in% c('both', 'ds', 'vs'))

axon_enriched_list <- as.list(axon_enriched$gene)

library(topGO)
library(ALL)
data(ALL)
data(geneList)
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)
sum(topDiffGenes(geneList))
sampleGOdata <- new("topGOdata", 
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)
sampleGOdata
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
resultKS.elim

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, 
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
allRes

pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]

gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)

plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
      elim = pValue.elim[sel.go],
      classic = pValue.classic[sel.go])
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo ='all')
