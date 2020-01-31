#Set anno to mouse
anno_version <- anno_mmus

#KW2 MB IP vs Input

name = 'mb_ip_vs_mb_ub_kw_batch_2'

metadata_filter = "region == 'mb' & batch_trap %in% c('b2') & outlier == 'n' & ip_pool == 'n'"

design <- ~ age + ovx + ip

min_counts = 10

reference_column <- 'ip'
reference_level <- 'input'

run_deseq()

#PK MB IP vs Input

name = 'mb_ip_vs_mb_ub_pk_batch_1'

metadata_filter = "region == 'mb' & batch_trap %in% c('b3') & outlier == 'n'"

design <- ~ age + ip

min_counts = 10

reference_column <- 'ip'
reference_level <- 'input'

run_deseq()

### **Midbrain IP combined vs Midbrain UB**

name = 'mb_ip_vs_mb_ub_combined'

metadata_filter = "region == 'mb' & batch_trap %in% c('b2', 'b3') & outlier == 'n' & ip_pool == 'n'"

design <- ~ batch_trap + age + ovx + ip

min_counts = 10

reference_column <- 'ip'
reference_level <- 'input'

run_deseq()

#Test agreement between KW and PK midbrain data
pk_mb_ip <- res_output_mb_ip_vs_mb_ub_pk_batch_1 %>% 
  filter(padj < 0.01 & log2FoldChange > 0) 

kw_mb_ip <- res_output_mb_ip_vs_mb_ub_kw_batch_2 %>% 
  filter(padj < 0.01 & log2FoldChange > 0) 

pk_mb_ip %>% 
  filter(gene %in% kw_mb_ip$gene) %>%
  nrow()

kw_mb_ip %>%
  filter(gene %in% pk_mb_ip$gene) %>%
  nrow()

### **Dorsal Striatum IP vs Dorsal Striatum UB**

name = 'ds_ip_vs_ds_ub'

metadata_filter = "region == 'ds' & outlier == 'n'"

design = ~ age + ip

min_counts = 10

reference_column <- 'ip'
reference_level <- 'input'

run_deseq()

### **Ventral Striatum IP vs Ventral Striatum UB**

name = 'vs_ip_vs_vs_ub'

metadata_filter = "region == 'vs' & outlier == 'n'"

design = ~ age + ip

min_counts = 10

reference_column <- 'ip'
reference_level <- 'input'

run_deseq()