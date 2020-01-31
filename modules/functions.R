run_blindedCounts <- function() {
  #subset sample metadata
  sample_metadata <- sample_metadata %>% filter(eval(parse(text=metadata_filter)))
  assign(paste('sample_metadata', name, sep = '_'), envir = .GlobalEnv, sample_metadata)
  
  ##prepare data for pca on all samples
  name = name
  design = ~ 1
  min_counts = 10
  
  #create list of kallisto input files
  files <- file.path(kallisto_path, sample_metadata$sample_code, "abundance.h5")
  names(files) <- sample_metadata$sample_name
  #import kallisto files
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
  #create deseq2 object
  dds <- DESeqDataSetFromTximport(txi, colData = sample_metadata, design = design)
  #filter low counts
  keep_feature <- rowMeans(counts(dds)) >= min_counts
  dds_removed <- dds[!keep_feature, ]
  dds <- dds[keep_feature, ]
  #generate dds and counts
  dds <- DESeq(dds, minReplicatesForReplace = replicate_replace)
  counts <- counts(dds, normalized = T) #DESeq normalises raw counts using median of ratios method, based on size factors
  vsd <- DESeq2::vst(dds, blind = T) #Blind for QC
  
  # assign(paste('dds_no_design', name, sep = '_'), envir = .GlobalEnv, dds)
  assign(paste('counts_blinded', name, sep = '_'), envir = .GlobalEnv, counts)
  assign(paste('vsd_blinded', name, sep = '_'), envir = .GlobalEnv, vsd)
  
}

#DESeq2 2020_01_28
run_deseq <- function() {
  #Is annotation human or mouse?
  anno <- anno_version
    #subset sample metadata
  sample_metadata <- sample_metadata %>% filter(eval(parse(text=metadata_filter)))
  assign(paste('sample_metadata', name, sep = '_'), envir = .GlobalEnv, sample_metadata)
  
  #create list of kallisto input files
  files <- file.path(kallisto_path, sample_metadata$sample_code, "abundance.h5")
  names(files) <- sample_metadata$sample_name
  
  #import kallisto files
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
  
  #create deseq2 object
  dds <- DESeqDataSetFromTximport(txi, colData = sample_metadata, design = design)
  
  #relevel contrast
  dds[[paste(reference_column)]] <- relevel(dds[[paste(reference_column)]], ref = reference_level)
  
  #filter low counts
  keep_feature <- rowMeans(counts(dds)) >= min_counts
  dds_removed <- dds[!keep_feature, ]
  dds <- dds[keep_feature, ]
  paste('Dimensions of sample data: ', dim(dds))
  
  #find ERCC spike-ins
  spikes <- rownames(dds)[grep("^ERCC", rownames(dds))]
  
  #calculate XYZ
  dds <- DESeq(dds, minReplicatesForReplace = replicate_replace)
  counts <- counts(dds, normalized = T)
  vsd <- DESeq2::vst(dds, blind = F)
  assign(paste('dds', name, sep = '_'), envir = .GlobalEnv, dds)
  assign(paste('dds_removed', name, sep = '_'), envir = .GlobalEnv, dds_removed)
  assign(paste('counts', name, sep = '_'), envir = .GlobalEnv, counts)
  assign(paste('vsd', name, sep = '_'), envir = .GlobalEnv, vsd)
  
  #summary(normalizationFactors(dds))
  
  #perform RUVseq without DESeq2 normalisation
  #eda_set <- newSeqExpressionSet(round(counts(dds, normalized = F)), phenoData = data.frame(sample_metadata, row.names=colnames(counts)))
  #EDASeq::plotRLE(eda_set, outline = F, ylim = c(-3, 3), col = colors()[sample_metadata$ip])
  
  #perform RUVseq with DESeq2 normalisation
  eda_set <- newSeqExpressionSet(round(counts(dds, normalized = T)), phenoData = data.frame(sample_metadata, row.names=colnames(counts)))
  EDASeq::plotRLE(eda_set, outline = F, ylim = c(-3, 3))
  
  # eda_set1 <- betweenLaneNormalization(eda_set, which="full")
  # plotRLE(eda_set1, outline = F, ylim = c(-3, 3), col = colors()[sample_metadata$ip])
  
  # eda_set1 <- RUVg(eda_set, spikes, k=1)
  # pData(eda_set1)
  # plotRLE(eda_set1, outline = F, ylim = c(-3, 3), col = colors()[sample_metadata$ip])
  
  #create contrast variable
  contrast_name <- tail(resultsNames(dds), n=1)
  contrast_split <- strsplit(contrast_name, '_')
  contrast <- c(contrast_split[[1]][1], contrast_split[[1]][2], contrast_split[[1]][4])
  
  #calculate result
  res <- results(dds, alpha = alpha, lfcThreshold = res_lfc, contrast = contrast, independentFiltering = T)
  assign(paste('res', name, sep = '_'), envir = .GlobalEnv, res)
  summary(res)
  
  #perform apeglm shrinkage
  res_shrunken <- lfcShrink(dds, coef = tail(resultsNames(dds), n=1), type = "apeglm", 
                            lfcThreshold = res_lfc, res = res, quiet = T)
  DESeq2::plotMA(res_shrunken, ylim=c(-5,5))
  assign(paste('res_shrunken', name, sep = '_'), envir = .GlobalEnv, res_shrunken)
  
  #add independent hypothesis weighting using basemean
  deRes <- as.data.frame(res)
  ihwRes <- ihw(pvalue ~ baseMean,  data = deRes, alpha = alpha)
  resIHW <- results(dds, filterFun=ihw, alpha = alpha, lfcThreshold = res_lfc, 
                    contrast = contrast)
  assign(paste('resIHW', name, sep = '_'), envir = .GlobalEnv, resIHW)
  
  res = resIHW
  
  #put results QC in (pvalue rejections, basemeans)
  
  #create results output object
  res_tb <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>%
    na.omit()
  
  sig <- res_tb %>%
    filter(padj < alpha & abs(log2FoldChange) > lfc)
  sig <- merge(as.data.frame(sig), anno, by.x = c("gene"), by.y = c("ensembl_gene_id"))
  sig <- merge(as.data.frame(sig), as.data.frame(counts(dds, normalized=TRUE)), 
               by.x = c("gene"), by.y = "row.names", sort=F)
  assign(paste('sig', name, sep = '_'), envir = .GlobalEnv, sig)
  
  res_tb <- merge(as.data.frame(res_tb), anno, by.x = c("gene"), by.y = c("ensembl_gene_id"))
  res_tb <- merge(as.data.frame(res_tb), as.data.frame(counts(dds, normalized=TRUE)), 
                  by.x = c("gene"), by.y = "row.names", sort=F)
  res_tb$padj[res_tb$padj == 0] <- .Machine$double.xmin #replacing absolute 0 values with smallest possible double for GSEA rank calculation
  assign(paste('res_output', name, sep = '_'), envir = .GlobalEnv, res_tb)
  
  res_tb$fcsign <- sign(res_tb$log2FoldChange)
  res_tb$logP <- -log10(res_tb$padj)
  res_tb$metric <- res_tb$logP/res_tb$fcsign
  res_tb <- res_tb %>% unique() %>% na.omit() %>% arrange(desc(metric)) %>% as_tibble
  rank <- res_tb[,c("external_gene_name", "metric")]
  assign(paste('rank', name, sep = '_'), envir = .GlobalEnv, rank)
  
  write.csv(res_tb, file = paste(file.path(output_dir), '/', 'results', '_', name, '.csv', sep = ''), row.names = F)
  write.table(rank, file = paste(file.path(output_dir), '/', 'rank', '_', name, '.rnk', sep = ''), sep = '\t', quote=F, row.names=F)
}

#PCA
make_pca <- function(name, keep_var, batch_correct) {
  if (batch_correct == T) {
    rv <- rowVars(assay(eval(as.name(paste('vsd_unblinded', name, batch, 'corrected', sep = '_'))))) #Calculate row variance
    rv_order <- order(rv, decreasing = TRUE)
    keep <- head(rv_order, max(1, nrow(assay(eval(as.name(paste('vsd_unblinded', name, batch, 'corrected', sep = '_')))))*(keep_var)))
    matrix_high_var <- assay(eval(as.name(paste('vsd_unblinded', name, batch, 'corrected', sep = '_'))))[keep, ] %>% #Select and transpose top n vsd genes
      t()
    pca <- prcomp(matrix_high_var, scale=T) #Calculate PCs
    var_explained <- pca$sdev^2/sum(pca$sdev^2) #Calculate PC variance
    assign(paste("pca_unblinded", name, batch, 'corrected', sep = '_'), envir = .GlobalEnv, pca)
    assign(paste("var_explained_unblinded", name, batch, 'corrected', sep = '_'), envir = .GlobalEnv, var_explained)
  }
  else {
    rv <- rowVars(assay(eval(as.name(paste("vsd_blinded", name, sep = "_"))))) #Calculate row variance
    rv_order <- order(rv, decreasing = TRUE)
    keep <- head(rv_order, max(1, nrow(assay(eval(as.name(paste('vsd_blinded', name, sep = '_')))))*(keep_var)))
    matrix_high_var <- assay(eval(as.name(paste("vsd_blinded", name, sep = "_"))))[keep, ] %>% #Select and transpose top n vsd genes
      t()
    pca <- prcomp(matrix_high_var, scale=T) #Calculate PCs
    var_explained <- pca$sdev^2/sum(pca$sdev^2) #Calculate PC variance
    assign(paste("pca_blinded", name, sep = '_'), envir = .GlobalEnv, pca)
    assign(paste("var_explained_blinded", name, sep = '_'), envir = .GlobalEnv, var_explained)
  }
}

plot_pca <- function(pca_df, group, label, title) {
  if (label == T) {
    group <- enquo(group)
    var_explained <- pca_df$sdev^2/sum(pca_df$sdev^2) #Calculate PC variance
    pca_df$x %>% 
      as.data.frame() %>%
      rownames_to_column("sample_name") %>%
      inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
      ggplot(aes(x=PC1,y=PC2, label = sample_name, color = !! group)) + 
      geom_point(aes(color = !! group), size = 4) +
      labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
           y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
      geom_label_repel(box.padding = 0.5) +
      theme_bw() +
      scale_color_lancet() + 
      scale_fill_lancet() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
      ggtitle(title) +
      theme(legend.position = "none")
  } else {
    group <- enquo(group)
    var_explained <- pca_df$sdev^2/sum(pca_df$sdev^2) #Calculate PC variance
    pca_df$x %>% 
      as.data.frame() %>%
      rownames_to_column("sample_name") %>%
      inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
      ggplot(aes(x=PC1,y=PC2, color = !! group)) + 
      geom_point(aes(color = !! group), size = 4) +
      labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
           y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
      theme_bw() +
      scale_color_lancet() + 
      scale_fill_lancet() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
      ggtitle(title) +
      theme(legend.position = "none")
  }
}

plot_scree <- function(pca_df){
  var_explained <- (pca_df$sdev^2/sum(pca_df$sdev^2))*100
  barplot(var_explained, 
          cex.names=1, 
          xlab=paste("Principal component (PC), 1-", 
                     length(pca$sdev)), 
          ylab="Proportion of variation (%)", 
          main="Scree plot", 
          ylim=c(0, 100))
}

plot_density <- function(vsd) {
  counts <- assay(vsd) %>% 
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    gather(key = sample_name, value = count, -gene) %>%
    inner_join(eval(as.name(paste("sample_metadata", name, sep = "_"))), by = c("sample_name" = "sample_name"))
  
  if (facet == T) {
    ggplot(counts, aes(count, fill = sample_name, colour = density_outlier)) +
      geom_line(stat = "density", adjust = 0.75) +
      theme_bw() +
      scale_color_lancet() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ylab("Density of Read Counts") +
      xlab("Variance Stabilised Read Count") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      facet_grid(rows = vars(eval(as.name(facet_factor))))
  } else {
    ggplot(counts, aes(count, fill = sample_name, colour = density_outlier)) +
      geom_line(stat = "density", adjust = 0.75) +
      theme_bw() +
      scale_color_lancet() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      ylab("Density of Read Counts") +
      xlab("Variance Stabilised Read Count") 
  }
}

plot_dendro <- function(vsd) {
  dendro <- as.dendrogram(hclust(d = dist(scale(t(assay(vsd))))))
  ggdendrogram(dendro, rotate = F)
}

#Produce unblinded VST counts
run_unblindedVst <- function() {
  #subset sample metadata
  sample_metadata <- sample_metadata %>% filter(eval(parse(text=metadata_filter)))
  assign(paste('sample_metadata', name, sep = '_'), envir = .GlobalEnv, sample_metadata)
  
  #create list of kallisto input files
  files <- file.path(kallisto_path, sample_metadata$sample_code, "abundance.h5")
  names(files) <- sample_metadata$sample_name
  #import kallisto files
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
  #create deseq2 object
  dds <- DESeqDataSetFromTximport(txi, colData = sample_metadata, design = design)
  #filter low counts
  keep_feature <- rowMeans(counts(dds)) >= min_counts
  dds_removed <- dds[!keep_feature, ]
  dds <- dds[keep_feature, ]
  #generate dds and counts
  dds <- DESeq(dds, minReplicatesForReplace = replicate_replace)
  counts <- counts(dds, normalized = T)
  vsd <- DESeq2::vst(dds, blind = F)
  
  assign(paste('vsd_unblinded', name, batch, 'corrected', sep = '_'), envir = .GlobalEnv, vsd)
  
}
