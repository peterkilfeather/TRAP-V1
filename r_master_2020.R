#Set up directories
base_dir <- getwd()
output_dir <- file.path(base_dir, 'output',paste('analysis', Sys.Date(), format(Sys.time(), "%H:%M"), sep = '_'))
dir.create(output_dir, recursive = T)
knitr::opts_knit$set(root.dir = output_dir)

kallisto_path <- file.path("/zfs/analysis/pk_trap/kallisto_201911_PK1_KW2/") #Ensembl v98 + human chromosome 4. Aligned November 2019 (See Github)

#Packages to load
library(DESeq2)
library(tximport)
# library(GenomicFeatures)
# library(GenomicAlignments)
# library(GenomicRanges)
# library(Rsamtools)
# library(Rsubread)
# library(refGenome)
# library(biomaRt)
# library(pheatmap)
# library(RColorBrewer)
library(RUVSeq)
library(IHW)
library(ggsci)
# library(WGCNA)
library(rhdf5)
library(naniar)
library(ggdendro)
# library(pcaExplorer)
# library(readxl)
library(ggrepel)
library(ggbeeswarm)
library(cowplot)
library(biomaRt)
library(AnnotationDbi)
# library(ensembldb)
# library(AnnotationHub)
library(Rhdf5lib)
# library(tidyverse)
# library(affy)
library(ggpubr)
library(gprofiler2)
# library(metap)
library(viridis)

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)
library(stringr)
library(forcats)
library(purrr)


filter <- dplyr::filter
select <- dplyr::select

#Functions to load
# source("objects/SampleNetwork_1.06.r")
# source("/zfs/analysis/r_functions/unblinded_vst_20191031.R")
# source("/zfs/analysis/r_functions/density_plot_20190905.R")
# source("/zfs/analysis/r_functions/dopa_counts_20191030.R")
# source("/zfs/analysis/r_functions/alignment_plot_20190906.R")
source(file.path(base_dir, "modules/functions.R"))

#Load annotation data
# txdb <- makeTxDbFromGFF('/zfs/analysis/pk_trap/kallisto_201911_PK1_KW2/reference/Mus_musculus.GRCm38.Homo_sapiens_GRCh38_combined_ens98_20191115.gtf', organism = "Mus musculus")
# saveDb(txdb, file = file.path(base_dir, "objects/Mus_musculus.GRCm38.Homo_sapiens_GRCh38_combined_ens98_20191115.sqlite"))
txdb <- loadDb(file.path(base_dir, "objects/Mus_musculus.GRCm38.Homo_sapiens_GRCh38_combined_ens98_20191115.sqlite"))
tx2gene <- AnnotationDbi::select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")
paste("The tx2gene table contains ", length(tx2gene$TXNAME), " entries.", sep = "")

anno_mmus <- readRDS(file.path(base_dir, "objects/anno_mmus.RDS")) %>%
  select(ensembl_gene_id, external_gene_name) %>%
  distinct()
anno_hsap <- readRDS(file.path(base_dir, "objects/anno_hsap.RDS")) %>%
  select(ensembl_gene_id, external_gene_name) %>%
  distinct()

ensembl_mmus <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# ensembl_hsap <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
# anno_hsap = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
# filters = 'ensembl_transcript_id', values = tx_hsap$TXNAME, mart = ensembl_hsap)

tx_mmus <- AnnotationDbi::select(txdb, keys(txdb, keytype = "TXNAME"), "TXNAME", "TXNAME") %>%
  as_tibble() %>%
  filter(str_detect(TXNAME, "^ENSMUST"))

tx_hsap <- AnnotationDbi::select(txdb, keys(txdb, keytype = "TXNAME"), "TXNAME", "TXNAME") %>%
  as_tibble() %>%
  filter(str_detect(TXNAME, "^ENST"))

#Close any open hd5 files
h5closeAll()

#Load sample metadata and QC statistics
p_pseudoaligned <- read_delim(file.path('/zfs/analysis/pk_trap/kallisto_201911_PK1_KW2/pseudoalignment_stats.txt'), delim = '\t', 
                              col_names = c("sample_code", "percentage_pseudoaligned"), col_types = "cd")

sample_metadata <- read_delim(file.path(base_dir, 'metadata/sample_metadata_20200128.csv'), delim = ",") %>%
  mutate(sample_code = as.character(sample_code)) %>%
  inner_join(p_pseudoaligned, by = c("sample_code" = "sample_code"))%>%
  mutate(align_rank = rank(percentage_pseudoaligned) / length(percentage_pseudoaligned))

# sample_metadata_kw_batch_2 <- sample_metadata %>%
#   as_tibble() %>%
#   filter(ip_pool == 'n' & batch_trap == 'b2')
# 
# sample_metadata_pk_batch_1_striatum <- sample_metadata %>% 
#   as_tibble() %>%
#   filter(ip_pool == 'n' & batch_trap == 'b3' & region != 'mb')

#Set up global parameters
skip_pcaExplorer = T
alpha = 0.01
res_lfc = 0 #lfc threshold for DESeq2 results functions
lfc = 0 #lfc threshold for post-results filtering genes
replicate_replace = Inf
dpi = 300
width = 7.25 #Dimensions set to span entire width of A4 publication, 4:3 ratio
height = 5.4375
units = "in"

#Run coverage ratio analysis
source(file.path(base_dir, "modules/coverage_ratio.R"))

#Quality control
source(file.path(base_dir, "modules/qc.R"))

#IP vs Input
source(file.path(base_dir, "modules/ip_vs_input.R"))

#Region Comparison
source(file.path(base_dir, "modules/region_diff_exp.R"))

