#Midbrain WGCNA

name = 'mb_18m_wt_vs_ovx_KW2'

metadata_filter = "region == 'mb' & ip == 'ip' & age == 'old' & batch_trap == 'b2' & ip_pool == 'n'"
metadata_filter = "region == 'mb' & ip == 'ip' & age == 'old' & batch_trap == 'b2' & ip_pool == 'n' & outlier == 'n'"

design = ~ sex + collection_day + ovx

min_counts = 100 #Rationale?

reference_column <- 'ovx'
reference_level <- 'wt'

run_deseq()

plotPCA(vsd_mb_18m_wt_vs_ovx_KW2, intgroup = 'ovx')

#Generate blinded counts
run_blindedCounts()
make_pca(name, keep_var = 1, batch_correct = F)
plot_pca(pca_blinded_mb_18m_wt_vs_ovx_KW2, group = ovx, label = T, title = 'Pre-WGCNA')

#Try and adjust for collection day
batch <- "collection_day"
min_counts <- 100
run_unblindedVst()

make_pca(name, keep_var = 1, batch_correct = T)
plot_pca(pca_unblinded_mb_18m_wt_vs_ovx_KW2_collection_day_corrected, group = ovx, label = T, title = 'Pre-WGCNA')

#Checking for good genes
gsg = goodSamplesGenes(assay(dds_mb_18m_wt_vs_ovx_KW2), verbose = 5);
gsg$allOK

vsd_mb_18m_wt_vs_ovx_KW2_corrected <- vst(dds_mb_18m_wt_vs_ovx_KW2, blind = F)
mat <- assay(vsd_mb_18m_wt_vs_ovx_KW2_corrected)
mat <- limma::removeBatchEffect(mat, vsd_mb_18m_wt_vs_ovx_KW2_corrected$collection_day)
assay(vsd_mb_18m_wt_vs_ovx_KW2_corrected) <- mat

plotPCA(vsd_mb_18m_wt_vs_ovx_KW2_corrected, intgroup = 'ovx')

datExprT <- assay(vsd_unblinded_mb_18m_wt_vs_ovx_KW2_collection_day_corrected)
datExprT <- as.data.frame(datExprT) %>%
  rownames_to_column(var = "gene")

all(sample_metadata_mb_18m_wt_vs_ovx_KW2$sample_name == colnames(datExprT[2:17]))

sample_metadata_mb_18m_wt_vs_ovx_KW2 <- as.data.frame(sample_metadata_mb_18m_wt_vs_ovx_KW2[-c(3,4)])

SampleNetwork(
  datExprT=datExprT,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=1,
  indices1=list(c(2:17)),
  sampleinfo1=sample_metadata_mb_18m_wt_vs_ovx_KW2,
  subgroup1=13,
  samplelabels1=1,
  grouplabels1=4,
  fitmodels1=TRUE,
  whichmodel1="multivariate",
  whichfit1="pc1",
  btrait1=c(3,6,13),
  trait1=c(13),
  asfactors1=c(3,6,13),
  projectname1="mb_18m_wt_vs_ovx",
  cexlabels1=0.7,
  normalize1=FALSE, #Don't quantile normalise DESeq2-normalised RNA-Seq data
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

#IP4F18mPD Z.K < -2

input = as.data.frame(t(datExprT[, -c(1)]))
names(input)<- datExprT$gene
dim(input)

traits <- sample_metadata_mb_18m_wt_vs_ovx_KW2[c(13)]
traits$ovx <- -1 * (as.numeric(factor(traits$ovx)) - 2)
rownames(traits) <- rownames(input)

#WGCNA Analysis
powers = c(c(1:11), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(input, corFnc = "bicor", corOptions = list(maxPOutliers=0.1), networkType = "signed hybrid", powerVector = powers, verbose = 5)

# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power (Signed)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

power <- 6

##AUTOMATIC METHOD!

datExpr <- input
datTraits <- traits

net = blockwiseModules(datExpr, 
                       power = power,
                       networkType = "signed hybrid",
                       TOMType = "signed",
                       corType = "bicor",
                       maxPOutliers = 0.1,
                       pearsonFallback = "individual",
                       maxBlockSize = 20000,
                       minModuleSize = 30, 
                       deepSplit = 2,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 5
)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]], 
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

moduleTraitCorInteresting <- as.data.frame(moduleTraitCor) %>%
  rownames_to_column("module") %>%
  filter(ovx > 0.5 | ovx  < -0.5) %>%
  column_to_rownames("module") %>%
  as.matrix()

moduleTraitPvalueInteresting <- as.data.frame(moduleTraitPvalue) %>%
  rownames_to_column("module") %>%
  filter(ovx < 0.1 & module %in% rownames(moduleTraitCorInteresting)) %>%
  column_to_rownames("module") %>%
  as.matrix()

# moduleTraitInteresting <- moduleTraitCorInteresting %>%
#   inner_join(y=moduleTraitPvalueInteresting, by = "module", suffix = c(".cor", ".pval")) %>%
#   column_to_rownames("module") %>%
#   as.matrix()


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCorInteresting, 2), "\n(",
                   signif(moduleTraitPvalueInteresting, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorInteresting)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCorInteresting,
               xLabels = names(datTraits),
               yLabels = rownames(moduleTraitCorInteresting),
               ySymbols = rownames(moduleTraitCorInteresting),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable ovx containing the ovx column of datTrait
ovx = as.data.frame(datTraits$ovx);
names(ovx) = "ovx"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, ovx, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(ovx), sep="");
names(GSPvalue) = paste("p.GS.", names(ovx), sep="");

module = "orange"
column = match(module, modNames)
moduleGenes = moduleColors==module # GET THIS WORKING TO HAVE CLEAR GENE NAMES, NOT ENSEMBL

module_mm_sig <-  tibble(module_membership = abs(geneModuleMembership[moduleGenes, column]), trait_significance = abs(geneTraitSignificance[moduleGenes, 1]))

module_plot <- ggplot(module_mm_sig, aes(x = module_membership, y = trait_significance)) + 
  geom_point() + 
  theme_classic() + 
  xlab('Module Membership') +
  ylab('Correlation with OVX Status') +
  expand_limits(y = 0) +
  scale_color_npg() +
  scale_fill_npg() +
  labs(fill = 'Region',
       title = 'Correlation of Module Membership with OVX Status',
       caption = "Genes situated in the top right of the graph are most likely to be of interest.") +
  theme(plot.title = element_text(hjust = 0, face = "bold"), plot.caption = element_text(hjust = 0, face = "bold"))

module_plot + 
  stat_cor(method = "pearson")

module_vsd <- input[, moduleGenes] %>%
  scale() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  as_tibble() %>%
  gather(key = gene, value = vsd_counts, -sample) %>%
  inner_join(sample_metadata_mb_18m_wt_vs_ovx_KW2, by = c("sample" = "sample_name")) %>%
  select(sample, gene, ovx, vsd_counts) %>%
  arrange(desc(ovx)) %>%
  mutate(ovx = factor(ovx, levels = c("wt", "ovx")))

wt <- module_vsd %>%
  filter(ovx == "wt") %>% 
  ggplot(mapping = aes(x = sample, y = gene, fill = vsd_counts)) +
    geom_tile() +
    # xlab(label = "Sample") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # facet_grid(cols = vars(ovx), scales = "free_x", space = "free_x", labeller = labeller(ovx = c(ovx = "OVX", wt = "WT"))) +
    scale_fill_gradient2(low = "green", mid = "black", high = "red") +
    ggtitle("Wild-type") +
    # theme(strip.placement = "outside", plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(),
          # strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
    # scale_y_discrete(limits = levels(as.factor(module_vsd$gene))) + 
    theme_cowplot() +
    theme(axis.text.y = element_blank()) + 
    theme(axis.ticks.y = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.4)) +
    theme(legend.position = "None") +
    theme(axis.line = element_blank()) + 
    theme(axis.title.x = element_blank()) +
    ylab("Module gene") +
    theme(axis.title.y = element_text(face = "bold", angle = 90))

ovx <- module_vsd %>%
  filter(ovx == "ovx") %>% 
  ggplot(mapping = aes(x = sample, y = gene, fill = vsd_counts)) +
  geom_tile() +
  # xlab(label = "Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # facet_grid(cols = vars(ovx), scales = "free_x", space = "free_x", labeller = labeller(ovx = c(ovx = "OVX", wt = "WT"))) +
  scale_fill_gradient2(low = "green", mid = "black", high = "red") +
  ggtitle("OVX") +
  # theme(strip.placement = "outside", plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(),
  # strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  # scale_y_discrete(limits = levels(as.factor(module_vsd$gene))) + 
  theme_cowplot() +
  theme(axis.text.y = element_blank()) + 
  theme(axis.ticks.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4)) +
  theme(legend.position = "None") +
  theme(axis.line = element_blank()) + 
  theme(axis.title.x = element_blank())
  
plot_row <- plot_grid(wt, ovx,
                      align = "hv",
                      rel_widths = c(0.55, 0.45))

# now add the title
title <- ggdraw() + 
  draw_label(
    "RNA processing module expression change in 18 month old OVX mice",
    fontface = 'bold',
    x = 0,
    hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot

x_lab <- ggdraw() + 
  draw_label(
    "Sample",
    fontface = 'bold',
    x = 0.5,
    hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(title, plot_row, x_lab,
          ncol = 1,
          rel_heights = c(0.1, 1, 0.1),
          align = "hv")

ggsave(path = output_dir, filename = "module_heatmap.png", width = width, height = height, dpi = dpi) 


# dev.off()
# which.module="orange";
# plotMat(t(scale(datExpr[,moduleGenes]) ),nrgcols=30,rlabels=T,
#         clabels=T,rcols=which.module,
#         title=which.module )


module_of_interest_genes <- module_vsd %>%
  inner_join(anno_mmus, by = c("gene"  = "ensembl_gene_id")) %>%
  select(external_gene_name) %>%
  arrange() %>%
  distinct() %>%
  pull()


background_genes <- rownames(assay(vsd_blinded_mb_18m_wt_vs_ovx_KW2)) %>%
  enframe(name = NULL) %>%
  rename("ensembl_gene_id" = "value") %>%
  inner_join(anno_mmus, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  select(external_gene_name) %>%
  pull()


gost(query = module_of_interest_genes, 
     organism = "mmusculus", 
     ordered_query = FALSE, 
     multi_query = FALSE, 
     significant = FALSE, 
     exclude_iea = T, 
     measure_underrepresentation = FALSE, 
     evcodes = FALSE, 
     user_threshold = 0.05, 
     correction_method = "gSCS", 
     domain_scope = "custom_annotated", 
     custom_bg = background_genes, 
     numeric_ns = "", 
     sources = "GO:BP", 
     as_short_link = T)






# # for the second (blue) module we use
# which.module="blue";
# plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
#         clabels=T,rcols=which.module,
#         title=which.module )
# which.module="brown";
# plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
#         clabels=T,rcols=which.module,
#         title=which.module )
# #write.csv(file = "ovx_module_1.csv", names(datExpr)[moduleColors=="brown"])


####Manual?

sampleTree = hclust(dist(input), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

n_genes = ncol(input)
n_samples = nrow(input)

traits <- sample_metadata_mb_18m_wt_vs_ovx[c(3, 6, 13)]
traits$collection_day <- as.numeric(traits$collection_day)-1
traits$sex <- as.numeric(traits$sex)-1
traits$ovx <- as.numeric(traits$ovx)-1
rownames(traits) <- rownames(input)

# Re-cluster samples
sampleTree2 = hclust(dist(input), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traits, signed = T);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traits),
                    main = "Sample dendrogram and trait heatmap")

adjacency = adjacency(input, power = power, type = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05, pearsonFallback = "individual"));

TOM = TOMsimilarity(adjacency, TOMType = "signed", verbose = 3);

dissTOM = 1-TOM

nSelect = 900
# For reproducibility, we set the random seed
set.seed(10);
select = sample(n_genes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, main = "Topological Overlap Matrix of 18 month IP samples")

geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

CUT_HEIGHT <- 0.995

abline(h = CUT_HEIGHT, col = "red");
clust = cutreeStatic(sampleTree, cutHeight = CUT_HEIGHT, minSize = 10)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize,
                            method = "hybrid");
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(input, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(input, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

#Calculate module eigengenes and correlate eigengenes to trait
MEs0 = moduleEigengenes(input, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, n_samples);

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

moduleTraitCor1 <- as.data.frame(moduleTraitCor) %>%
  rownames_to_column("module") %>%
  filter(ovx > 0.5 | ovx < -0.5)

moduleTraitPvalue1 <- as.data.frame(moduleTraitPvalue) %>%
  rownames_to_column("module") %>%
  filter(ovx < 0.1)

# Define variable X containing the X column of traits
ovx = as.data.frame(traits$ovx);
names(ovx) = "ovx"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(input, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), n_samples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(input, ovx, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), n_samples));
names(geneTraitSignificance) = paste("GS.", names(ovx), sep="");
names(GSPvalue) = paste("p.GS.", names(ovx), sep="");

module = "coral"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for disease",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")


which.module="coral";
plotMat(input[,colorh1==which.module], nrgcols=30,
        clabels=rownames(input),rcols=which.module)

datME=moduleEigengenes(heatmap_INPUT,colorh1)$eigengenes

sizeGrWindow(8,7);
which.module="brown"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
plotMat(t(scale(heatmap_INPUT[,colorh1==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module, cex.main=1, clabels=rownames(heatmap_INPUT))
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression")


