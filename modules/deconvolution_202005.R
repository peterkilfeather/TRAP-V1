# Plot DA Markers
genes <- c('ENSMUSG00000000214', #Th 
           'ENSMUSG00000021609', #Slc6a3
           'ENSMUSG00000025094', #Vmat2
           'ENSMUSG00000026826') #Nurr1 

name = 'axon_ip_vs_ub'
metadata_filter = "batch_trap %in% c('b3') & region %in% c('ds', 'vs') & outlier == 'n'"
run_blindedCounts()

all_counts <- as.data.frame(counts_blinded_axon_ip_vs_ub + 1)
anno <- anno_mmus

counts_to_plot <- all_counts %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble %>%
  filter(gene %in% genes) %>% 
  gather(sample,counts,-gene) %>%
  inner_join(anno, by = c("gene" = "ensembl_gene_id")) %>%
  inner_join(sample_metadata, by = c("sample" = "sample_name")) %>%
  mutate(ip = recode(ip, ip = "IP", input = "Input"))

title <- 'Dopaminergic markers in axon TRAP samples'
caption <- 'Markers of dopaminergic neurons are enriched in axonal samples.'
caption <- paste0(strwrap(caption, 120), sep="", collapse="\n")

ip_ub_multi_boxplot <- ggplot(counts_to_plot, aes(x = external_gene_name, y = counts, fill = ip)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge()) +
  theme_bw() +
  xlab('Gene') +
  ylab('log2 Counts') +
  scale_color_lancet() +
  scale_fill_lancet() +
  labs(fill = 'Sample Type',
       title = title, 
       caption = caption) +
  theme(plot.caption = element_text(hjust = 0, face = "bold"), legend.position = "right")

ip_ub_multi_boxplot




striatal_markers <- read_delim("external_datasets/striatal_markers.csv", delim = ",") %>%
  rename("cell_type" = `Cell Type `) %>%
  inner_join(anno_mmus, by = c('gene' = 'external_gene_name'))

res_output <- res_output_ds_ip_vs_ds_ub

res_output %>%
  inner_join(striatal_markers, by = c('gene' = 'ensembl_gene_id')) %>%
  mutate(signif = ifelse(padj < 0.01, TRUE, FALSE)) %>%
  arrange(desc(baseMean)) %>%
  ggplot(aes(x = cell_type, y = log2FoldChange)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

cell_types = marker_res %>% pull(cell_type) %>% unique()


fisher <- function(a,b,c,d){
  data <- matrix(c(a,b,c,d),ncol=2)
  c(p = fisher.test(data, alternative = "two.sided")$p.value,
    OR = fisher.test(data, alternative = "two.sided")$estimate)
}

res_output %>%
  inner_join(striatal_markers, by = c('gene' = 'ensembl_gene_id')) %>%
  mutate(signif = ifelse(padj < 0.01, TRUE, FALSE)) %>%
  arrange(desc(baseMean)) %>%
  select(gene, baseMean, log2FoldChange, padj, external_gene_name, cell_type, signif) %>%
  group_by(cell_type) %>%
  summarise(de_go = sum(signif),
            go = n(), 
            mean_lfc = mean(log2FoldChange)) %>%
  mutate(de_no_go = dim(res_output %>% filter(padj < 0.01))[1] - de_go,
         no_de_go = go - de_go,
         no_de_no_go = dim(res_output)[1] - dim(res_output %>% filter(padj < 0.01))[1] - no_de_go) %>% 
  rowwise() %>%
  mutate(p=fisher(de_go, no_de_go, de_no_go, no_de_no_go)[[1]],
         OR=fisher(de_go, no_de_go, de_no_go, no_de_no_go)[[2]]) #Two-sided probability

phyper(27, 
       46, 
       dim(res_output_ds_ip_vs_ds_ub)[1] - 46, 
       dim(res_output_ds_ip_vs_ds_ub %>% filter(padj < 0.01))[1], 
       lower.tail = FALSE) #probability of getting 28 or more genes of GO term

sum(dhyper(28:46, 
       46, 
       dim(res_output_ds_ip_vs_ds_ub)[1] - 46, 
       dim(res_output_ds_ip_vs_ds_ub %>% filter(padj < 0.01))[1])) #probability of getting 28 or more genes of GO term

phyper(1, 
       21, 
       dim(res_output_ds_ip_vs_ds_ub)[1] - 21, 
       dim(res_output_ds_ip_vs_ds_ub %>% filter(padj < 0.01))[1]) #probability of getting 1 or fewer genes of GO term

sum(dhyper(0:1, 
       21, 
       dim(res_output_ds_ip_vs_ds_ub)[1] - 21, 
       dim(res_output_ds_ip_vs_ds_ub %>% filter(padj < 0.01))[1])) #probability of getting 1 or fewer genes of GO term
           
                  
anno_mmus_upper <- anno_mmus %>%
  mutate(external_gene_name_upper = toupper(external_gene_name))

cell_markers_mckenzie <- read_delim("external_datasets/cell_type_markers_mckenzie.csv", delim = ",") %>% 
  inner_join(anno_mmus_upper, by = c("gene" = "external_gene_name_upper")) %>%
  select(gene, cell_type, ensembl_gene_id)
  
res_output %>%
  inner_join(cell_markers_mckenzie, by = c('gene' = 'ensembl_gene_id')) %>%
  mutate(signif = ifelse(padj < 0.01, TRUE, FALSE)) %>%
  arrange(desc(baseMean)) %>%
  select(gene, baseMean, log2FoldChange, padj, external_gene_name, cell_type, signif) %>%
  group_by(cell_type) %>%
  summarise(de_go = sum(signif),
            go = n(), 
            mean_lfc = mean(log2FoldChange)) %>%
  mutate(de_no_go = dim(res_output %>% filter(padj < 0.01))[1] - de_go,
         no_de_go = go - de_go,
         no_de_no_go = dim(res_output)[1] - dim(res_output %>% filter(padj < 0.01))[1] - no_de_go) %>% 
  rowwise() %>%
  mutate(p=fisher(de_go, no_de_go, de_no_go, no_de_no_go)[[1]],
         OR=fisher(de_go, no_de_go, de_no_go, no_de_no_go)[[2]]) #Two-sided probability

res_output %>%
  inner_join(cell_markers_mckenzie, by = c('gene' = 'ensembl_gene_id')) %>%
  mutate(signif = ifelse(padj < 0.01, TRUE, FALSE)) %>%
  arrange(desc(baseMean)) %>%
  ggplot(aes(x = cell_type, y = log2FoldChange)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_cowplot()




# 
# name = 'pk_batch_1_ip_str_no_outliers'
# metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region %in% c('ds', 'vs') & outlier == 'n'"
# run_blindedCounts()
# 
# str_ip_percents <- counts_blinded_pk_batch_1_ip_str_no_outliers %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   as_tibble() %>%
#   pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
#   group_by(gene) %>%
#   summarise(mean_count = mean(count)) %>%
#   mutate(percent_total = mean_count*100/sum(mean_count)) %>%
#   arrange(desc(percent_total)) %>%
#   inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id")) %>% 
#   mutate(rank = 1:n())
# 
# name = 'pk_batch_2_ub_str'
# metadata_filter = "batch_trap %in% c('b3') & ip == 'input' & region %in% c('ds', 'vs')"
# run_blindedCounts()
# 
# str_ub_percents <- counts_blinded_pk_batch_2_ub_str %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   as_tibble() %>%
#   pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
#   group_by(gene) %>%
#   summarise(mean_count = mean(count)) %>%
#   mutate(percent_total = mean_count*100/sum(mean_count)) %>%
#   arrange(desc(percent_total)) %>%
#   inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id")) %>% 
#   mutate(rank = 1:n())

str_fold_change <- counts_blinded_axon_ip_vs_ub %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  inner_join(sample_metadata, by = c('sample' = 'sample_name')) %>%
  group_by(gene, ip) %>%
  summarise(mean_count = mean(count)) %>%
  select(gene, ip, mean_count) %>%
  pivot_wider(names_from = ip, values_from = mean_count) %>%
  mutate(fold_change = ip/input) %>%
  inner_join(anno_mmus, by = c('gene' = 'ensembl_gene_id')) %>%
  ungroup() %>%
  mutate(input_percent = input*100/sum(input), 
         ip_percent = ip*100/sum(ip),
         percent_fold_change = ip_percent/input_percent) %>%
  arrange(desc(fold_change)) %>%
  mutate(fc_rank = 1:n()) 


str_fold_change <- counts_blinded_axon_ip_vs_ub %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  inner_join(sample_metadata, by = c('sample' = 'sample_name')) %>%
  group_by(gene, ip) %>%
  summarise(mean_count = mean(count)) %>%
  select(gene, ip, mean_count) %>%
  pivot_wider(names_from = ip, values_from = mean_count) %>%
  mutate(fold_change = ip/input) %>%
  inner_join(anno_mmus, by = c('gene' = 'ensembl_gene_id')) %>%
  ungroup() %>%
  mutate(input_percent = input*100/sum(input), 
         ip_percent = ip*100/sum(ip),
         percent_fold_change = ip_percent/input_percent) %>%
  arrange(desc(fold_change)) %>%
  mutate(fc_rank = 1:n()) 



name = 'mb_ip_vs_ub'
metadata_filter = "batch_trap %in% c('b2') & region %in% c('mb') & outlier == 'n'"
run_blindedCounts()


mb_fold_change <- counts_blinded_mb_ip_vs_ub %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  inner_join(sample_metadata, by = c('sample' = 'sample_name')) %>%
  group_by(gene, ip) %>%
  summarise(mean_count = mean(count)) %>%
  select(gene, ip, mean_count) %>%
  pivot_wider(names_from = ip, values_from = mean_count) %>%
  mutate(fold_change = ip/input) %>%
  inner_join(anno_mmus, by = c('gene' = 'ensembl_gene_id')) %>%
  ungroup() %>%
  mutate(input_percent = input*100/sum(input), 
         ip_percent = ip*100/sum(ip),
         percent_fold_change = ip_percent/input_percent) %>%
  arrange(desc(fold_change)) %>%
  mutate(fc_rank = 1:n()) 

mb_fold_change %>% 
  ggplot(aes(fold_change)) +
  geom_line(stat = "density", adjust = 0.75) +
  coord_cartesian(xlim=c(0, 10)) +
  geom_vline(xintercept = 1)

changes_str <- str_ip_percents %>%
  inner_join(str_ub_percents, by = c('gene' = 'gene'), suffix = c(".ip", ".ub"))

mb_ip_percents <- counts_blinded_pk_batch_1_ip_mb_round_4 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  group_by(gene) %>%
  summarise(mean_count = mean(count)) %>%
  mutate(percent_total = mean_count*100/sum(mean_count)) %>%
  arrange(desc(percent_total)) %>%
  inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id")) %>% 
  mutate(rank = 1:n())

name = 'kw_batch_2_ub_mb'
metadata_filter = "batch_trap %in% c('b2') & ip == 'input' & region == 'mb'"
run_blindedCounts()

mb_ub_percents <- counts_blinded_kw_batch_2_ub_mb %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  group_by(gene) %>%
  summarise(mean_count = mean(count)) %>%
  mutate(percent_total = mean_count*100/sum(mean_count)) %>%
  arrange(desc(percent_total)) %>%
  inner_join(anno_mmus, by = c("gene" = "ensembl_gene_id")) %>% 
  mutate(rank = 1:n())

changes_mb <- mb_ip_percents %>%
  inner_join(mb_ub_percents, by = c('gene' = 'gene'), suffix = c(".ip", ".ub"))

temp <- changes_mb %>%
  mutate(rank_change = rank.ip-rank.ub) %>%
  arrange(rank_change)

temp <- changes_str %>%
  mutate(rank_change = rank.ip-rank.ub) %>%
  arrange(rank_change)




# %>%
#   slice(3:100) %>%
#   ggplot(aes(x = fct_reorder(gene, desc(percent_total)), y = percent_total)) +
#   geom_point() +
#   theme_cowplot()




for (cell in cell_types){
  # Initialize variables
  m <- dim(marker_res %>% filter(cell_type == cell))[1]       # Genes IN GO term
  n <- dim(res_output_ds_ip_vs_ds_ub)[1] - m # Genes NOT IN GO term
  k <- dim(res_output_ds_ip_vs_ds_ub %>% filter(padj < 0.01))[1]       # Gene hits, that is, differentially expressed
  x <- c(0:m)  # Genes both IN GO term and differentially expressed 'hits'
  
  # Use the dhyper built-in function for hypergeometric density
  probabilities <- dhyper(x, m, n, k, log = FALSE)
  
  # Calculate the one-sided p-value for ? or more genes both DE and IN GO term.
  pvalue <- sum(probabilities[dim(marker_res %>% filter(cell_type == cell & signif == 'y'))[1]:dim(marker_res %>% filter(cell_type == cell))[1]])
  pvalue
  
  

  # Bar plot
  data <- data.frame( x = x, y = probabilities )
  print(ggplot(data, aes(x=factor(x), y=y)) +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"),
          axis.title.x=element_text(margin=margin(20,0,0,0)),
          axis.title.y=element_text(margin=margin(0,20,0,0))
    ) +
    geom_bar(stat="identity", fill=ifelse(data$x < dim(marker_res %>% filter(cell_type == cell & signif == 'y'))[1],
                                          rgb(52, 73, 94, maxColorValue=255),
                                          rgb(231, 76, 60, maxColorValue=255)),
             colour="black") +
    labs(x = "DE genes IN GO term", y = "Probability") +
      ggtitle(cell) + 
    theme_cowplot()
  )
  
  
}













library(SRAdb)
# srafile = getSRAdbFile()
sra_dbname <- file.path("/zfs/analysis/SRAmetadb.sqlite")
sra_con <- dbConnect( dbDriver("SQLite"), sra_dbname)
ascpCMD <- 'ascp -QT -l 300m -i /home/peter/.aspera/connect/bin/ascp'

sraConvert("SRX4111655", sra_con)
sraConvert(c('SRX4111659'), sra_con= sra_con)

listSRAfile(in_acc = c("SRP001007"), sra_con = sra_con, fileType = 'fastq', srcType='fasp')










genes <- c('ENSMUSG00000030088', #aldh1l1
           'ENSMUSG00000020932', #gfap
           'ENSMUSG00000033208', #s100b
           'ENSMUSG00000005089', #Slc1a2/GLT-1
           'ENSMUSG00000005360') #Slc1a3/GLAST

counts_to_plot <- all_counts %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble %>%
  filter(gene %in% genes) %>% 
  gather(sample,counts,-gene) %>%
  inner_join(anno, by = c("gene" = "ensembl_gene_id")) %>%
  inner_join(sample_metadata, by = c("sample" = "sample_name")) %>%
  mutate(ip = recode(ip, ip = "IP", input = "Input"))

title <- 'Astrocyte markers in axon TRAP samples'
caption <- 'Markers of astroctytes are enriched in axonal samples, except S100b.'
caption <- paste0(strwrap(caption, 120), sep="", collapse="\n")

ip_ub_multi_boxplot <- ggplot(counts_to_plot, aes(x = external_gene_name, y = counts, fill = ip)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge()) +
  theme_bw() +
  xlab('Gene') +
  ylab('log2 Counts') +
  scale_color_npg() +
  scale_fill_npg() +
  labs(fill = 'Sample Type',
       title = title, 
       caption = caption) +
  theme(plot.caption = element_text(hjust = 0, face = "bold"), legend.position = "right")

ip_ub_multi_boxplot
