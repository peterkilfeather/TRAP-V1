# **Dorsal Striatum IP vs Midbrain IP**

name = 'ds_ip_vs_mb_ip'

metadata_filter = "ip == 'ip' & outlier == 'n' & batch_trap %in% c('b3') & region %in% c('ds', 'mb')"

design <- ~ collection_day + age + region

min_counts = 10

reference_column <- 'region'
reference_level <- 'mb'

run_deseq()

# **Ventral Striatum IP vs Midbrain IP**

name = 'vs_ip_vs_mb_ip'

metadata_filter = "ip == 'ip' & outlier == 'n' & batch_trap %in% c('b3') & region %in% c('vs', 'mb')"

design <- ~ collection_day + age + region

min_counts = 10

reference_column <- 'region'
reference_level <- 'mb'

run_deseq()

# **Identifying region specific genes - October 2019**

enrichment_cutoff = 3 #1 = 2-fold

#Collect genes enriched in ip
ds_up_ip <- res_output_ds_ip_vs_ds_ub %>%
  filter(padj < 0.01 & log2FoldChange > 0)
vs_up_ip <- res_output_vs_ip_vs_vs_ub %>%
  filter(padj < 0.01 & log2FoldChange > 0)
axon_up_ip <- ds_up_ip %>%
  filter(gene %in% vs_up_ip$gene)
mb_up_ip <- res_output_mb_ip_vs_mb_ub_combined %>%
  filter(padj < 0.01 & log2FoldChange > 0)

#Collect genes depleted in ip
ds_down_ip <- res_output_ds_ip_vs_ds_ub %>%
  filter(padj < 0.01 & log2FoldChange < 0)
vs_down_ip <- res_output_vs_ip_vs_vs_ub %>%
  filter(padj < 0.01 & log2FoldChange < 0)
axon_down_ip <- ds_down_ip %>%
  filter(gene %in% vs_down_ip$gene)
mb_down_ip <- res_output_mb_ip_vs_mb_ub_combined %>%
  filter(padj < 0.01 & log2FoldChange < 0)

#Collect genes regionally enriched (ds/vs/mb)

ds_vs_mb_enriched <- res_output_ds_ip_vs_mb_ip %>%
  filter(gene %in% ds_up_ip$gene & log2FoldChange > enrichment_cutoff & padj < 0.01)
vs_vs_mb_enriched <- res_output_vs_ip_vs_mb_ip %>%
  filter(gene %in% vs_up_ip$gene & log2FoldChange > enrichment_cutoff & padj < 0.01)
axon_vs_mb_enriched <- ds_vs_mb_enriched %>%
  filter(gene %in% vs_vs_mb_enriched$gene)
mb_vs_ds_enriched <- res_output_ds_ip_vs_mb_ip %>%
  filter(gene %in% mb_up_ip$gene & log2FoldChange < -enrichment_cutoff & padj < 0.01)
mb_vs_vs_enriched <- res_output_vs_ip_vs_mb_ip %>%
  filter(gene %in% mb_up_ip$gene & log2FoldChange < -enrichment_cutoff & padj < 0.01)
mb_vs_axon_enriched <- mb_vs_ds_enriched %>% filter(gene %in% mb_vs_vs_enriched$gene)

# write.table(ds_vs_mb_enriched, file = paste(file.path(output_dir), '/', 'ds_vs_mb_enriched.txt', sep = ''))
# write.table(vs_vs_mb_enriched, file = paste(file.path(output_dir), '/', 'vs_vs_mb_enriched.txt', sep = ''))
# write.table(axon_vs_mb_enriched, file = paste(file.path(output_dir), '/', 'axon_vs_mb_enriched.txt', sep = ''))
# write.table(mb_vs_ds_enriched, file = paste(file.path(output_dir), '/', 'mb_vs_ds_enriched.txt', sep = ''))
# write.table(mb_vs_vs_enriched, file = paste(file.path(output_dir), '/', 'mb_vs_vs_enriched.txt', sep = ''))
# write.table(mb_vs_axon_enriched, file = paste(file.path(output_dir), '/', 'mb_vs_axon_enriched.txt', sep = ''))


# 8-fold Striatum vs Midbrain, 2-fold Dorsal vs Ventral for Richard and Natalie, 2019_11_28

name = 'ds_ip_vs_vs_ip'

metadata_filter = "axon == 'axon' & outlier == 'n' & ip == 'ip'"

design = ~ age + region

min_counts = 10

reference_column <- 'region'
reference_level <- 'ds'

run_deseq()


ds_vs_vs_vs_mb_enriched <- res_output_ds_ip_vs_vs_ip  %>%
  filter(gene %in% ds_up_ip$gene & log2FoldChange > 1 & padj < 0.01) %>%
  filter(gene %in% ds_vs_mb_enriched$gene)

goi <- ds_vs_vs_vs_mb_enriched %>%
  select(external_gene_name, gene, baseMean, log2FoldChange, padj)

goi_ds_vs_mb <- ds_vs_mb_enriched %>%
  filter(gene %in% goi$gene) %>%
  select(external_gene_name, gene, baseMean, log2FoldChange, padj)

greater_in_ds_than_vs_and_mb <- goi %>%
  inner_join(goi_ds_vs_mb, by = c("gene" = "gene")) %>%
  select("gene_name" = external_gene_name.x, 
         "ensembl_id" = gene,
         "base_mean_ds_vs_vs" = baseMean.x, 
         "log2_foldchange_ds_vs_vs" = log2FoldChange.x,
         "adjusted_p_ds_vs_vs" = padj.x,
         "base_mean_ds_vs_mb" = baseMean.y, 
         "log2_foldchange_ds_vs_mb" = log2FoldChange.y,
         "adjusted_p_ds_vs_mb" = padj.y) %>%
  write_csv(file.path(output_dir, "DS_selective.csv"))

## Volcano plot
res_output_ds_ip_vs_mb_ip_plot <- res_output_ds_ip_vs_mb_ip %>%
  arrange(desc(padj)) %>%
  mutate(goi = ifelse(gene %in% greater_in_ds_than_vs_and_mb$ensembl_id, 'GOI', 'Background'))


ggplot(res_output_ds_ip_vs_mb_ip_plot, 
       aes(x = log2FoldChange, y =-log10(padj), colour = goi)) +
  geom_point(alpha = 0.3) + 
  geom_vline(xintercept = 1, linetype="dotted") +
  geom_label_repel(aes(label = ifelse(gene %in% greater_in_ds_than_vs_and_mb$ensembl_id, as.character(external_gene_name), '')), 
                   box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') + 
  ggtitle("Regional Differences in Gene Expression") + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  theme_bw() +
  scale_color_lancet()  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# 
# ## Volcano plot
# res_output_axon_ip_vs_mb_ip_plot_2 <- res_output_axon_ip_vs_mb_ip %>%
#   arrange(desc(padj)) %>%
#   mutate(Axonal = ifelse(padj < 0.01 & log2FoldChange > 3, 'Axonal', 'Somal')) %>%
#   mutate(Axonal = factor(Axonal, levels = c("Somal", "Axonal"))) %>%
#   mutate(plot_label = ifelse(padj < 0.000001 & log2FoldChange > 2, 'y', 'n')) %>%
#   mutate(goi = ifelse(gene %in% greater_in_ds_than_vs_and_mb$ensembl_id, 'GOI', 'Background'))
# 
# 
# ggplot(res_output_axon_ip_vs_mb_ip_plot_2, 
#        aes(x = log2FoldChange, y =-log10(padj), colour = Axonal)) +
#   geom_point(alpha = 0.3) + 
#   geom_vline(xintercept = 3, linetype="dotted") +
#   geom_label_repel(aes(label = ifelse(gene %in% greater_in_ds_than_vs_and_mb$ensembl_id, as.character(external_gene_name), '')), 
#                    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') + 
#   ggtitle("Dorsal Striatum vs Midbrain") + 
#   xlab("log2 fold change") + 
#   ylab("-log10 adjusted p-value") +
#   #scale_y_continuous(limits = c(0,50)) +
#   theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25))) +
#   theme_bw() +
#   scale_color_lancet()  +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())
# 

res_output_ds_ip_vs_vs_ip_plot <- res_output_ds_ip_vs_vs_ip %>%
  mutate(Dorsal = ifelse(padj < 0.01 & log2FoldChange > 1, 'Dorsal', 'Background')) %>%
  mutate(goi = ifelse(gene %in% greater_in_ds_than_vs_and_mb$ensembl_id, 'GOI', 'Background')) %>%  
  mutate(goi = factor(goi, levels= c("Background", "GOI")))

ggplot(res_output_ds_ip_vs_vs_ip_plot, 
       aes(x = log2FoldChange, y =-log10(padj), colour = Dorsal)) +
  geom_point(alpha = 0.3) + 
  # geom_vline(xintercept = 1, linetype="dotted") +
  geom_label_repel(aes(label = ifelse(gene %in% greater_in_ds_than_vs_and_mb$ensembl_id, as.character(external_gene_name), '')), 
                   box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') + 
  ggtitle("Dorsal Striatum vs Ventral Striatum") + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  theme_bw() +
  scale_color_lancet()  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(path = output_dir, filename = "dorsal_vs_ventral_goi.png") 


# 8-fold Dorsal vs Ventral, regardless of MB differential for Richard and Natalie, 2019_11_28

goi_1lfc <- res_output_ds_ip_vs_vs_ip %>%
  filter(gene %in% ds_up_ip$gene & log2FoldChange > 1 & padj < 0.01) %>%
  select(external_gene_name, gene, baseMean, log2FoldChange, padj)
write_csv(goi_1lfc, file.path(output_dir, "ds_vs_vs_2fold.csv"))

goi_2lfc <- res_output_ds_ip_vs_vs_ip %>%
  filter(gene %in% ds_up_ip$gene & log2FoldChange > 2 & padj < 0.01) %>%
  select(external_gene_name, gene, baseMean, log2FoldChange, padj)
write_csv(goi_2lfc, file.path(output_dir, "ds_vs_vs_4fold.csv"))

goi_3lfc <- res_output_ds_ip_vs_vs_ip %>%
  filter(gene %in% ds_up_ip$gene & log2FoldChange > 3 & padj < 0.01) %>%
  select(external_gene_name, gene, baseMean, log2FoldChange, padj)
write_csv(goi_3lfc, file.path(output_dir, "ds_vs_vs_8fold.csv"))

res_output_ds_ip_vs_vs_ip %>%
  filter(gene %in% ds_up_ip$gene) %>%
  mutate(Dorsal = ifelse(gene %in% goi_3lfc$gene, '8-fold', 
                         ifelse(gene %in% goi_2lfc$gene, '4-fold',
                                ifelse(gene %in% goi_1lfc$gene, '2-fold', 'Background')))) %>%
  ggplot(aes(x = log2FoldChange, y =-log10(padj), colour = Dorsal)) +
  geom_point(alpha = 0.3) + 
  # geom_vline(xintercept = 1, linetype="dotted") +
  geom_label_repel(aes(label = ifelse(gene %in% goi_3lfc$gene, as.character(external_gene_name),
                                      ifelse(gene %in% goi_2lfc$gene, as.character(external_gene_name), ''))), 
                   box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') + 
  ggtitle("Dorsal Striatum vs Ventral Striatum") + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  theme_cowplot() +
  scale_color_lancet()  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  coord_cartesian(xlim = c(-2.5, 5), 
                  ylim = c(0, 100)) +
  theme(axis.title = element_text(face = "bold"))

ggsave(file.path(output_dir,  "ds_vs_vs_volcano_RWM_NCR_20201128.png"), width = width, height = height, dpi = dpi, units = units) 

