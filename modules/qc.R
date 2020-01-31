#PCA all IP samples
name = 'pk_batch_1_all_ip'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)
sample_metadata <- sample_metadata %>% 
  mutate(round_1_pca_outlier = ifelse(sample_name == "IPMB18MC", 'y', 'n'))

plot_pca(pca_blinded_pk_batch_1_all_ip, group = region, label = F, title = 'All samples')
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_1_region.png", width = width, height = height, dpi = dpi, units = units) 
plot_pca(pca_blinded_pk_batch_1_all_ip, group = round_1_pca_outlier, label = F, title = 'All samples')
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_1_outiers.png", width = width, height = height, dpi = dpi, units = units) 

#MB IP: IPMB18MC Outlier Evidence
name = 'pk_batch_1_ip_all_mb'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region == 'mb'"
run_blindedCounts()

#Density plot of old MB IP samples, highlighting IPMB18MC
assay(vsd_blinded_pk_batch_1_ip_all_mb) %>% 
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample_name, value = count, -gene) %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>% 
  filter(age == "old") %>%
  ggplot(aes(count, fill = sample_name, colour = round_1_pca_outlier)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_bw() +
  scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Density of Read Counts") +
  xlab("Variance Stabilised Read Count") +
  facet_grid(rows = vars(sample_name)) + 
  theme(legend.position = "none") + 
  theme(strip.text.y = element_text(angle = 360)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_1_IPMB18MC_density.png", width = width, height = height, dpi = dpi, units = units)

#PCA after IPMB18MC removal
name = 'pk_batch_1_ip_round_2'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & round_1_pca_outlier == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)

sample_metadata <- sample_metadata %>% 
  mutate(round_2_pca_outlier = ifelse(sample_name %in% c("IPDS3MB", "IPDS3ME", "IPDS18MC"), 'y', 'n'))

plot_pca(pca_blinded_pk_batch_1_ip_round_2, group = region, label = F, title = 'IPMB18MC Removed')
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_2_region.png", width = width, height = height, dpi = dpi, units = units) 

plot_pca(pca_blinded_pk_batch_1_ip_round_2, group = round_2_pca_outlier, label = F, title = 'IPMB18MC Removed')
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_2_outliers.png", width = width, height = height, dpi = dpi, units = units) 

#DS IP: IPDS3MB, IPDS3ME, IPDS18MC Outlier Evidence
name = 'pk_batch_1_ip_ds'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region == 'ds'"
run_blindedCounts()

sample_metadata <- sample_metadata %>% 
  mutate(round_2_density_outlier = ifelse(sample_name == "IPDS3MB", 'y', 'n'))
sample_metadata <- sample_metadata %>% 
  mutate(round_2_dopa_outlier = ifelse(sample_name %in% c("IPDS3ME", "IPDS18MC"), 'y', 'n'))

#Density plot of young DS IP, highlighting IPDS3MB
assay(vsd_blinded_pk_batch_1_ip_ds) %>% 
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample_name, value = count, -gene) %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>% 
  filter(age == "young") %>%
  ggplot(aes(count, fill = sample_name, colour = round_2_density_outlier)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_bw() +
  scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Density of Read Counts") +
  xlab("Variance Stabilised Read Count") +
  facet_grid(rows = vars(sample_name)) + 
  theme(legend.position = "none") + 
  theme(strip.text.y = element_text(angle = 360)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(path = output_dir, filename="pca_blinded_pk_batch_1_all_ip_round_2_IPDS3MB_density.png", width = width/1.5, height = height, dpi = dpi, units = units)

#Dopa counts highlighting IPDS3ME and IPDS18MC
assay(vsd_blinded_pk_batch_1_ip_ds) %>% 
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample_name, value = count, -gene) %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>% 
  filter(gene == 'ENSMUSG00000000214' | gene == 'ENSMUSG00000021609') %>%
  inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id")) %>% 
  ggplot(counts, mapping = aes(x = reorder(sample_name, count), fill = round_2_dopa_outlier)) + 
  geom_bar(width = 0.8, position = "dodge", aes(weight = count)) +
  theme_bw() +
  scale_color_lancet() +
  scale_fill_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) +
  ylab("Variance Stabilised Read Count") + 
  facet_grid(rows = vars(external_gene_name)) + 
  theme(legend.position = "none") + 
  theme(strip.text.y = element_text(angle = 360)) 

ggsave(path = output_dir, filename="pca_blinded_pk_batch_1_all_ip_round_2_IPDS3ME_IPDS18MC_dopa.png", width = width/1.5, height = height, dpi = dpi, units = units)

#PCA after IPMB18MC, IPDS3MB, IPDS3ME, IPDS18MC removal
name = 'pk_batch_1_ip_round_3'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & round_1_pca_outlier == 'n' & round_2_pca_outlier == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)

sample_metadata <- sample_metadata %>% 
  mutate(round_3_pca_outlier = ifelse(sample_name %in% c("IPVS3MC"), 'y', 'n'))

plot_pca(pca_blinded_pk_batch_1_ip_round_3, group = region, label = F, title = 'MB and DS Outliers Removed')
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_3_region.png", width = width, height = height, dpi = dpi, units = units) 

plot_pca(pca_blinded_pk_batch_1_ip_round_3, group = round_3_pca_outlier, label = F, title = 'MB and DS Outliers Removed')
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_3_outliers.png", width = width, height = height, dpi = dpi, units = units) 

#VS IP: IPVS3MC Outlier Evidence
name = 'pk_batch_1_ip_vs'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region == 'vs'"
run_blindedCounts()

#Density plot of young VS IP samples, showing no difference
assay(vsd_blinded_pk_batch_1_ip_vs) %>% 
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample_name, value = count, -gene) %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>% 
  filter(age == "young") %>%
  ggplot(aes(count, fill = sample_name, colour = round_3_pca_outlier)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_bw() +
  scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Density of Read Counts") +
  xlab("Variance Stabilised Read Count") +
  facet_grid(rows = vars(sample_name)) + 
  theme(legend.position = "none") + 
  theme(strip.text.y = element_text(angle = 360)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(path = output_dir, filename="pca_blinded_pk_batch_1_all_ip_round_3_IPVS3M_density.png", width = width/1.5, height = height, dpi = dpi, units = units)

#Dopa counts highlighting IPVS3MC
assay(vsd_blinded_pk_batch_1_ip_vs) %>% 
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample_name, value = count, -gene) %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>% 
  filter(gene == 'ENSMUSG00000000214' | gene == 'ENSMUSG00000021609') %>%
  inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id")) %>% 
  ggplot(counts, mapping = aes(x = reorder(sample_name, count), fill = round_3_pca_outlier)) + 
  geom_bar(width = 0.8, position = "dodge", aes(weight = count)) +
  theme_bw() +
  scale_color_lancet() +
  scale_fill_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) +
  ylab("Variance Stabilised Read Count") + 
  facet_grid(rows = vars(external_gene_name)) + 
  theme(legend.position = "none") + 
  theme(strip.text.y = element_text(angle = 360))

ggsave(path = output_dir, filename="pca_blinded_pk_batch_1_all_ip_round_3_IPVS3MC_dopa.png", width = width/1.5, height = height, dpi = dpi, units = units)

sample_filter <- c("IPVS18MD", "IPVS3MC")
assay(vsd_blinded_pk_batch_1_ip_vs) %>%
  cor() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  as_tibble() %>%
  filter(!sample %in% sample_filter) %>%
  replace_with_na_if(.predicate = is.numeric, 
                     condition = ~.x > 0.999) %>%
  pivot_longer(-sample, names_to = "correlated_sample", values_to = "corr") %>%
  drop_na(corr) %>%
  ggplot(aes(correlated_sample, corr)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_lancet() + 
  scale_fill_lancet() + 
  geom_quasirandom() +
  # geom_jitter(width = 0.2) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(path = output_dir, filename="corr_plot_pk_batch_1_vs_ip_round_3.png", width = width, height = height, dpi = dpi, units = units)

#PCA after IPMB18MC, IPDS3MB, IPDS3ME, IPDS18MC, IPVS3MC removal
name = 'pk_batch_1_ip_round_4'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & round_1_pca_outlier == 'n' & round_2_pca_outlier == 'n' & round_3_pca_outlier == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)

plot_pca(pca_blinded_pk_batch_1_ip_round_4, group = region, label = F, title = 'MB, DS and VS Outliers Removed')
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_4_region.png", width = width, height = height, dpi = dpi, units = units) 

#Demonstrate how choice of number of variables can influence clustering
dir.create(file.path(output_dir, "pca_movie"))
n <- 1
name = 'pk_batch_1_all_ip'

for (i in seq(0.002, 0.2, 0.002)) {
  make_pca(name, i, batch_correct = F)
  plot_pca(pca_blinded_pk_batch_1_all_ip, group = region, label = F, title = paste("PCA using top ", i*100, "% variable features", sep = ""))
  ggsave(path = file.path(output_dir, "pca_movie"), filename = sprintf("region%03d.png", n), width = width, height = height, dpi = dpi, units = units) 
  plot_pca(pca_blinded_pk_batch_1_all_ip, group = suspected_outlier, label = F, title = paste("PCA using top ", i*100, "% variable features", sep = ""))
  ggsave(path = file.path(output_dir, "pca_movie"), filename = sprintf("outlier%03d.png", n), width = width, height = height, dpi = dpi, units = units) 
  n <- n + 1
}

assay(eval(as.name(paste('vsd_blinded', name, sep = '_')))) %>%
  as.data.frame() %>%
  rownames_to_column("gene")%>%
  mutate(variance = rowVars(assay(eval(as.name(paste('vsd_blinded', name, sep = '_')))))) %>%
  mutate(mean = rowMeans(assay(eval(as.name(paste('vsd_blinded', name, sep = '_')))))) %>%
  mutate(include = ifelse(variance > quantile(variance, 0.95,), 'y', 'n')) %>%
  select(gene, mean, variance, include) %>%
  arrange(desc(variance)) %>%
  mutate(order = 1:n()) %>%
  ggplot(aes(mean, variance, colour = include, alpha = 1/2)) +
  geom_point()  +
  theme_cowplot() +
  scale_color_lancet() +
  ylab("Variance") +
  xlab("Mean Expression") + 
  theme(legend.position = "none")
ggsave(path = output_dir, filename="pca_blinded_pk_batch_1_all_ip_round_1_meanVariance_top5percent.png", width = width, height = height, dpi = dpi) 

#MB PCA after outliers removed
name = 'pk_batch_1_ip_mb_round_4'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region == 'mb' & round_1_pca_outlier == 'n' & round_2_pca_outlier == 'n' & round_3_pca_outlier == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)
plot_pca(pca_blinded_pk_batch_1_ip_mb_round_4, group = age, label = F, title = 'MB')

#adjust for batch
design <- ~ collection_day + age
batch <- "collection_day"
min_counts <- 10
run_unblindedVst()

make_pca(name, keep_var = 1, batch_correct = T)

#Demonstrate pairs plot problem
pairs(pca_blinded_pk_batch_1_ip_mb_round_4$x[,1:10], col="black", main="Principal components analysis bi-plot\nPCs 1-10", pch=16)

#Plot PCs on single axis
var_explained0 <- pca_blinded_pk_batch_1_ip_mb_round_4$sdev^2/sum(pca_blinded_pk_batch_1_ip_mb_round_4$sdev^2) #Calculate PC variance
names(var_explained0) <- colnames(pca_blinded_pk_batch_1_ip_mb_round_4$x)
var_explained <- as.data.frame(var_explained0) %>% 
  rownames_to_column('PC') %>%
  mutate(var = paste0((round(var_explained0*100)), ' %', sep = '')) %>%
  select(PC, var)

pca_blinded_pk_batch_1_ip_mb_round_4$x %>% 
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
  pivot_longer(cols = PC1:PC10, names_to = 'PC', values_to = 'x') %>%
  select(sample_name, PC, x, age) %>% 
  inner_join(var_explained, by = c('PC' = 'PC')) %>% 
  mutate(x_label = paste0(PC, '\n', var, sep = '')) %>%
  mutate(dummy = factor(PC, levels= c(1, 3, 4, 5, 6, 7, 8, 9, 10, 2))) %>%
  ggplot(aes(x = x_label, y = x, colour = age)) + 
  geom_point() +
  labs(x = paste0("Principle Component"),
       y = paste0("Component Coordinate")) +
  theme_bw() +
  scale_color_lancet() + 
  scale_fill_lancet() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  ggtitle("MB Principle Components") 
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_4_mb_single_axis.png", width = width, height = height, dpi = dpi) 

#DS PCA after outliers removed
name = 'pk_batch_1_ip_ds_round_4'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region == 'ds' & round_1_pca_outlier == 'n' & round_2_pca_outlier == 'n' & round_3_pca_outlier == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)
plot_pca(pca_blinded_pk_batch_1_ip_ds_round_4, group = age, label = F, title = 'DS')

# #adjust for batch
# design <- ~ collection_day + age
# batch <- "collection_day"
# min_counts <- 10
# run_unblindedVst()

# make_pca(name, keep_var = 1, batch_correct = T)

#Plot PCs on single axis
var_explained0 <- pca_blinded_pk_batch_1_ip_ds_round_4$sdev^2/sum(pca_blinded_pk_batch_1_ip_ds_round_4$sdev^2) #Calculate PC variance
names(var_explained0) <- colnames(pca_blinded_pk_batch_1_ip_ds_round_4$x)
var_explained <- as.data.frame(var_explained0) %>% 
  rownames_to_column('PC') %>%
  mutate(var = paste0((round(var_explained0*100)), ' %', sep = '')) %>%
  select(PC, var)

pca_blinded_pk_batch_1_ip_ds_round_4$x %>% 
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
  pivot_longer(cols = PC1:(length(var_explained0)+1), names_to = 'PC', values_to = 'x') %>%
  select(sample_name, PC, x, age) %>% 
  inner_join(var_explained, by = c('PC' = 'PC')) %>% 
  mutate(x_label = paste0(PC, '\n', var, sep = '')) %>%
  ggplot(aes(x = x_label, y = x, colour = age)) + 
  geom_point() +
  labs(x = paste0("Principle Component"),
       y = paste0("Component Coordinate")) +
  theme_bw() +
  scale_color_lancet() + 
  scale_fill_lancet() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  ggtitle("DS Principle Components") 
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_4_ds_single_axis.png", width = width, height = height, dpi = dpi) 

#VS PCA after outliers removed
name = 'pk_batch_1_ip_vs_round_4'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region == 'vs' & round_1_pca_outlier == 'n' & round_2_pca_outlier == 'n' & round_3_pca_outlier == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)
plot_pca(pca_blinded_pk_batch_1_ip_vs_round_4, group = age, label = T, title = 'VS Principle Components')

# #adjust for batch
# design <- ~ collection_day + age
# batch <- "collection_day"
# min_counts <- 10
# run_unblindedVst()

# make_pca(name, keep_var = 1, batch_correct = T)

#Plot PCs on single axis
var_explained0 <- pca_blinded_pk_batch_1_ip_vs_round_4$sdev^2/sum(pca_blinded_pk_batch_1_ip_vs_round_4$sdev^2) #Calculate PC variance
names(var_explained0) <- colnames(pca_blinded_pk_batch_1_ip_vs_round_4$x)
var_explained <- as.data.frame(var_explained0) %>% 
  rownames_to_column('PC') %>%
  mutate(var = paste0((round(var_explained0*100)), ' %', sep = '')) %>%
  select(PC, var)

pca_blinded_pk_batch_1_ip_vs_round_4$x %>% 
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
  pivot_longer(cols = PC1:(length(var_explained0)+1), names_to = 'PC', values_to = 'x') %>%
  select(sample_name, PC, x, age) %>% 
  inner_join(var_explained, by = c('PC' = 'PC')) %>% 
  mutate(x_label = paste0(PC, '\n', var, sep = '')) %>%
  ggplot(aes(x = x_label, y = x, colour = age)) + 
  geom_point() +
  labs(x = paste0("Principle Component"),
       y = paste0("Component Coordinate")) +
  theme_bw() +
  scale_color_lancet() + 
  scale_fill_lancet() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  ggtitle("VS Principle Components") 
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_4_vs_single_axis.png", width = width, height = height, dpi = dpi) 

#Remove IPVS18MD
sample_metadata <- sample_metadata %>% 
  mutate(round_4_pca_outlier = ifelse(sample_name %in% c("IPVS18MD"), 'y', 'n'))

name = 'pk_batch_1_ip_vs_round_5'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region == 'vs' & round_1_pca_outlier == 'n' & round_2_pca_outlier == 'n' & round_3_pca_outlier == 'n' & round_4_pca_outlier == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)
plot_pca(pca_blinded_pk_batch_1_ip_vs_round_5, group = age, label = T, title = 'VS Principle Components')

#Plot PCs on single axis
var_explained0 <- pca_blinded_pk_batch_1_ip_vs_round_5$sdev^2/sum(pca_blinded_pk_batch_1_ip_vs_round_5$sdev^2) #Calculate PC variance
names(var_explained0) <- colnames(pca_blinded_pk_batch_1_ip_vs_round_5$x)
var_explained <- as.data.frame(var_explained0) %>% 
  rownames_to_column('PC') %>%
  mutate(var = paste0((round(var_explained0*100)), ' %', sep = '')) %>%
  select(PC, var)

pca_blinded_pk_batch_1_ip_vs_round_5$x %>% 
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
  pivot_longer(cols = PC1:(length(var_explained0)+1), names_to = 'PC', values_to = 'x') %>%
  select(sample_name, PC, x, age) %>% 
  inner_join(var_explained, by = c('PC' = 'PC')) %>% 
  mutate(x_label = paste0(PC, '\n', var, sep = '')) %>%
  ggplot(aes(x = x_label, y = x, colour = age)) + 
  geom_point() +
  labs(x = paste0("Principle Component"),
       y = paste0("Component Coordinate")) +
  theme_bw() +
  scale_color_lancet() + 
  scale_fill_lancet() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  ggtitle("VS Principle Components") 
ggsave(path = output_dir, filename = "pca_blinded_pk_batch_1_all_ip_round_5_vs_single_axis.png", width = width, height = height, dpi = dpi) 

#UB PCA with demonstration of UBDS18MPOOLA-C dopa counts
name = 'pk_batch_1_all_ub'
metadata_filter = "batch_trap %in% c('b3') & ip == 'input'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)
sample_metadata <- sample_metadata %>% 
  mutate(round_1_ub_pca_outlier = ifelse(sample_name == "UBDS18MPOOLA-C", 'y', 'n')) 

plot_pca(pca_blinded_pk_batch_1_all_ub, group = region, label = F, title = 'All samples')

#Dopa counts highlighting UBDS18MPOOLA-C
assay(vsd_blinded_pk_batch_1_all_ub) %>% 
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample_name, value = count, -gene) %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>% 
  filter(region == 'ds') %>%
  filter(gene == 'ENSMUSG00000000214' | gene == 'ENSMUSG00000021609') %>%
  inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id")) %>% 
  ggplot(counts, mapping = aes(x = reorder(sample_name, count), fill = round_1_ub_pca_outlier)) + 
  geom_bar(width = 0.8, position = "dodge", aes(weight = count)) +
  theme_bw() +
  scale_color_lancet() +
  scale_fill_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) +
  ylab("Variance Stabilised Read Count") + 
  facet_grid(rows = vars(external_gene_name)) +
  theme(legend.position = "none") + 
  theme(strip.text.y = element_text(angle = 360)) 

#KW2 PCA
name = 'kw_batch_2_ip_mb'
metadata_filter = "batch_trap %in% c('b2') & ip == 'ip' & region == 'mb' & ip_pool == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)

sample_metadata <- sample_metadata %>% 
  mutate(round_1_KW2_pca_outlier = ifelse(sample_name == "IP3M3m", 'y', 'n'))

plot_pca(pca_blinded_kw_batch_2_ip_mb, group = age, label = T, title = 'KW2 MB')

plot_pca(pca_blinded_kw_batch_2_ip_mb, group = round_1_KW2_pca_outlier, label = T, title = 'KW2 MB')

#PCA after IPMB18MC removal
name = 'kw_batch_2_ip_mb_round_2'
metadata_filter = "batch_trap %in% c('b2') & ip == 'ip' & region == 'mb' & ip_pool == 'n' & round_1_KW2_pca_outlier == 'n'"
run_blindedCounts()

make_pca(name, keep_var = 1, batch_correct = F)

plot_pca(pca_blinded_kw_batch_2_ip_mb_round_2, group = region, label = T, title = 'IP3M3m Removed')
ggsave(path = output_dir, filename = "pca_blinded_KW_batch_2_ip_mb_round_2_region.png", width = width, height = height, dpi = dpi, units = units) 

# Produce final outlier list
sample_metadata <- sample_metadata %>%
  mutate(outlier = ifelse(round_1_pca_outlier == 'y', 'y', 
                          ifelse(round_2_pca_outlier == 'y', 'y', 
                                 ifelse(round_3_pca_outlier == 'y', 'y', 
                                        ifelse(round_4_pca_outlier == 'y', 'y', 
                                               ifelse(round_1_ub_pca_outlier == 'y', 'y',
                                                      ifelse(round_1_KW2_pca_outlier == 'y', 'y', 'n')))))))

sample_metadata %>% filter(outlier == 'y')

