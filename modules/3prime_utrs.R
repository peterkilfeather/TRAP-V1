
# **Generate 3 prime UTR fasta files for motif analysis**

# ids <- axon_vs_mb_enriched$gene
# axon_background_ids <- axon_up_ip$gene
# 
# axon_vs_mb_enriched_utr3 = getSequence(id = ids,
#                                        type = "ensembl_gene_id",
#                                        seqType = "3utr",
#                                        mart = ensembl_mmus) %>%
#   group_by(ensembl_gene_id) %>%
#   mutate(id = paste(ensembl_gene_id, row_number(), sep = '')) %>%
#   ungroup() %>%
#   select('3utr', id) %>%
#   rename('ensembl_gene_id' = id) %>%
#   as.data.frame()
# 
# axon_background_utr3 = getSequence(id = axon_background_ids,
#                                    type="ensembl_gene_id",
#                                    seqType="3utr",
#                                    mart=ensembl_mmus) %>%
#   group_by(ensembl_gene_id) %>%
#   mutate(id = paste(ensembl_gene_id, row_number(), sep = '')) %>%
#   ungroup() %>%
#   select('3utr', id) %>%
#   rename('ensembl_gene_id' = id) %>%
#   as.data.frame()
# 
# exportFASTA(axon_vs_mb_enriched_utr3, file.path(output_dir, 'axon_vs_mb_enriched_utr3.fa'))
# exportFASTA(axon_background_utr3, file.path(output_dir, 'axon_background_utr3.fa'))


#QAPA

salmon_path <- file.path("/zfs/analysis/pk_trap/salmon_qapa_20191123")

metadata <- filter(sample_metadata, ip == 'ip' & ip_pool == 'n' & outlier == 'n' & batch_trap == 'b3') %>% 
  select(c('sample_name', 'sample_code', 'age', 'region', 'axon'))
metadata <- metadata %>%
  mutate(path = file.path(salmon_path, metadata$sample_code, "quant.sf"))

sample_to_group <- metadata %>% select(sample_name, sample_code, age, region, axon)

#create list of input files
files <- file.path(salmon_path, metadata$sample_code, "quant.sf")
names(files) <- metadata$sample

##QAPA output
matrix0 <- read_tsv("/zfs/analysis/pk_trap/salmon_qapa_20191123/pau_results.txt")

matrix <- matrix0 %>%
  filter(Num_Events > 1) %>%
  replace(., is.na(.), 0)

apa_to_gene <- matrix %>%
  select(APA_ID, Gene_Name)

dim(matrix)
matrix$Gene %>% unique() %>% length() ##Number of genes
matrix$Transcript %>% unique() %>% length() ##Number of transcripts

#PCA from tpm
tpm <- select(matrix, APA_ID, contains("TPM")) %>%
  rename_at(vars(matches(".TPM")), ~ str_remove(., ".TPM")) %>%
  as_tibble()

tpm_metadata <- tpm %>%
  pivot_longer(-APA_ID, names_to = "sample_code", values_to = "TPM") %>%
  inner_join(metadata, by = c("sample_code" = "sample_code"))

tpm_subset <- tpm_metadata %>%
  select(APA_ID, sample_name, TPM) %>%
  pivot_wider(names_from = "sample_name", values_from = "TPM", 
              values_fn = list(TPM = mean))

tpm_matrix <- data.matrix(tpm_subset[-1])
rownames(tpm_matrix) <- tpm_subset$APA_ID
keep <- rowVars(tpm_matrix) > 0

matrix_high_var <- tpm_matrix[keep,] %>% #Select and transpose top n vsd genes
  t()
pca_pk_qapa_all_ip_tpm <- prcomp(matrix_high_var, scale=T) #Calculate PCs
var_explained <- pca_pk_qapa_all_ip_tpm$sdev^2/sum(pca_pk_qapa_all_ip_tpm$sdev^2) #Calculate PC variance

plot_pca(pca_pk_qapa_all_ip_tpm, group = region, label = F, title = 'All samples')

pc1_loadings0 <- pca_pk_qapa_all_ip_tpm$rotation %>%
  as.data.frame() %>% 
  rownames_to_column('APA_ID') %>%
  select(APA_ID, PC1) %>%
  mutate(PC1_abs = abs(PC1)) %>%
  arrange(desc(PC1_abs))

pc1_loadings <- filter(pc1_loadings0, PC1_abs > (quantile(pc1_loadings0$PC1_abs, 0.98))) %>%
  inner_join(apa_to_gene, by = c("APA_ID", "APA_ID")) %>%
  inner_join(anno, by = c("Gene_Name" = "Gene.name")) %>%
  select(APA_ID, Gene_Name, Gene.stable.ID, Gene.description, PC1)

ggplot(as.data.frame(pca_pk_qapa_all_ip_tpm$rotation[, 1:2]), aes(x = PC1, y = PC2)) + geom_point()

##TPM metrics
matrix_tpm_means <- tpm_metadata %>%
  inner_join(apa_to_gene, by = c("APA_ID" = "APA_ID")) %>%
  group_by(Gene_Name, sample_name) %>%
  summarise(tpm_sum_gene = sum(TPM)) %>% 
  group_by(Gene_Name) %>%
  summarise(tpm_mean_gene = mean(tpm_sum_gene)) 

genes_over_3_tpm <- matrix_tpm_means %>%
  filter(tpm_mean_gene > 3)
apa_over_3_tpm <- apa_to_gene %>%
  filter(Gene_Name %in% genes_over_3_tpm$Gene_Name)
apa_tpm_means <- tpm_metadata %>%
  group_by(APA_ID) %>%
  summarise(apa_tpm_mean = mean(TPM))

##PAU
#PCA from tpm
pau <- select(matrix, APA_ID, contains("PAU")) %>%
  rename_at(vars(matches(".PAU")), ~ str_remove(., ".PAU")) %>%
  as_tibble() %>%
  filter(APA_ID %in% apa_over_3_tpm$APA_ID) 

pau_renamed <- pau %>%
  pivot_longer(-APA_ID, names_to = "sample_code", values_to = "PAU") %>%
  inner_join(sample_metadata, by = c("sample_code" = "sample_code")) %>%
  filter(ip == 'ip') %>%
  select(APA_ID, sample_code, sample_name, PAU, age, axon) %>%
  group_by(axon, APA_ID) %>%
  summarise(mean = mean(PAU), median = median(PAU), min = min(PAU), max = max(PAU)) 
# %>% 
#   pivot_wider(id_cols = APA_ID, names_from = sample_name, values_from = PAU) 
  
  
  prox_pau_metadata <- pau %>%
  pivot_longer(-APA_ID, names_to = "sample_code", values_to = "PAU") %>%
  inner_join(metadata, by = c("sample_code" = "sample_code")) %>%
  filter(str_detect(APA_ID, '_P'))

prox_pau_subset <- prox_pau_metadata %>%
  select(APA_ID, sample_name, PAU) %>%
  pivot_wider(names_from = sample_name, values_from = PAU)

prox_pau_matrix <- data.matrix(prox_pau_subset[-1])
rownames(prox_pau_matrix) <- prox_pau_subset$APA_ID
keep <- rowVars(prox_pau_matrix) > 0

matrix_high_var <- prox_pau_matrix[keep,] %>% #Select and transpose top n vsd genes
  t()
pca_pk_qapa_all_ip_prox_pau <- prcomp(matrix_high_var, scale=T) #Calculate PCs

pca_df <- pca_pk_qapa_all_ip_prox_pau
var_explained <- pca_df$sdev^2/sum(pca_df$sdev^2) #Calculate PC variance
pca_df$x %>% 
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
  ggplot(aes(x = PC1, y = PC2, color = region)) + 
  geom_point(aes(color = axon), size = 4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme_cowplot() +
  scale_color_lancet() + 
  scale_fill_lancet() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  ggtitle('PCA using Proximal PolyA Usage Estimates') +
  labs(color = "Region")

ggsave(file.path(output_dir, "pca_ppau_region.png"), width = width, height = height, dpi = dpi) 


pc1_loadings0_prox_pau <- pca_pk_qapa_all_ip_prox_pau$rotation %>%
  as.data.frame() %>% 
  rownames_to_column('APA_ID') %>%
  select(APA_ID, PC1) %>%
  mutate(PC1_abs = abs(PC1)) %>%
  arrange(desc(PC1_abs))

# pc1_loadings_prox_pau_top50 <- pc1_loadings0_prox_pau %>%
#   top_n(50)
# 
# top50_prox_pau <- prox_pau_metadata %>%
#   filter(APA_ID %in% pc1_loadings_prox_pau_top50$APA_ID) %>%
#   group_by(APA_ID, axon) %>%
#   summarise(median_prox_pau = median(PAU)) %>%
#   pivot_wider(names_from = axon, values_from = median_prox_pau) %>%
#   mutate(diff = axon - mb) %>%
#   inner_join(pc1_loadings_prox_pau_top50, by = c("APA_ID" = "APA_ID"))  %>%
#   inner_join(apa_to_gene, by = c("APA_ID", "APA_ID")) %>%
#   inner_join(anno, by = c("Gene_Name" = "Gene.name")) %>%
#   inner_join(apa_tpm_means, by = c("APA_ID" = "APA_ID"))

prox_pau_region <- prox_pau_metadata %>%
  # filter(APA_ID %in% pc1_loadings_ppau_top50$APA_ID) %>%
  group_by(APA_ID, axon) %>%
  summarise(median_prox_pau = median(PAU)) %>%
  pivot_wider(names_from = axon, values_from = median_prox_pau) %>%
  mutate(diff = axon - mb)  %>% 
  inner_join(apa_to_gene, by = c("APA_ID", "APA_ID")) %>%
  inner_join(anno, by = c("Gene_Name" = "Gene.name")) %>%
  inner_join(apa_tpm_means, by = c("APA_ID" = "APA_ID")) %>%
  inner_join(matrix_tpm_means, by = c("Gene_Name" = "Gene_Name")) %>%
  inner_join(pc1_loadings0_prox_pau, by = c("APA_ID" = "APA_ID"))

ggplot(prox_pau_region, aes(x = diff, y = PC1)) + 
  geom_point() +
  theme_cowplot()


changes <- prox_pau_region %>%
  mutate(change = ifelse(diff > 10, 'axon_shortening', 
                         ifelse(diff < -10, 'axon_lengthening', 'no_change'))) %>%
  select(APA_ID, change) %>%
  group_by(change) %>%
  summarise(count = n())


prox_pau_region %>%
  mutate(change = ifelse(diff > 10, 'axon_shortening', 
                         ifelse(diff < -10, 'axon_lengthening', 'no_change'))) %>%
  ggplot(aes(x = diff)) +
  geom_area(stat = "density") +
  theme_cowplot()

1.96*sd(prox_pau_region$diff)

ggplot(changes, aes(change)) +
  geom_bar(aes(weight = count, fill = change)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_lancet() + 
  scale_fill_lancet() +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())


over_10 <- prox_pau_region %>%
  # mutate(shift = ifelse(APA_ID %in% axon_longer$APA_ID, 'longer', 'shorter')) %>%
  filter(abs(diff) > 10)

axon_longer <- over_10 %>%
  filter(diff < 0)
axon_longer_ip_enriched <- axon_longer %>%
  filter(Gene_Name %in% res_output_ds_ip_vs_ds_ub$Gene.name)

axon_shorter <- over_10 %>%
  filter(diff > 0)

ggsave(file.path(output_dir, "utr_lengthening_barchart.png"), width = width, height = height, dpi = dpi, units = units) 


ggplot(prox_pau_region, aes(x = diff, y = PC1, colour = PC1_abs, alpha = PC1_abs)) + 
  geom_point() +
  geom_vline(xintercept = 10, linetype = "dotted") +
  geom_vline(xintercept = -10, linetype = "dotted") +
  theme_cowplot() +
  scale_color_continuous(low = "#56B1F7", high = "#132B43") +
  scale_alpha_continuous(range = c(0.2, 1)) +
  # theme(panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank()) +
  xlab("\u0394 Proximal PolyA Usage (%)") +
  ylab("PC1 Loading") +
  ggtitle("3\' UTR Switching in Dopaminergic Axons") + 
  # labs(color='Axonal 3\'\nUTR Length') + 
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-60, 60)) +
  theme(axis.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "tandem_utr_switching.png"), width = width, height = height, dpi = dpi, units = units) 




