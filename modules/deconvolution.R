library(debCAM)

name = 'pk_batch_1_ip_str_no_outliers'
metadata_filter = "batch_trap %in% c('b3') & ip == 'ip' & region %in% c('ds', 'vs') & outlier == 'n'"
run_blindedCounts()

counts_blinded_pk_batch_1_ip_str_no_outliers

data <- counts_blinded_pk_batch_1_ip_ds_round_4

set.seed(81)

gene_means <- tibble(gene=names(rowMeans(counts_blinded_pk_batch_1_ip_ds_round_4)), 
                     mean_count = rowMeans(counts_blinded_pk_batch_1_ip_ds_round_4)) %>%
  arrange(desc(mean_count)) %>%
  mutate(percentile = rank(mean_count)/length(mean_count) * 100) 


quantile(gene_means$mean_count, probs = seq(0, 1, 0.05))

rCAM <- CAM(data, K = 2:7, 
            corner.strategy = 1,
            dim.rdc = 7, #There are 6 principal components generated in prcomp
            thres.low = 0.05, 
            thres.high = 0.95,
            cluster.method = "K-Means", #alternative is apcluster
            cluster.num = 100, #number of clusters to form. Must be >> k. default = 50
            lof.thres = 0.02, #local outlier removal using 'lofactor'. MG
            quick.select = 20,
            sample.weight = NULL,
            generalNMF = TRUE, #Default = FALSE
            cores = 8)

plot(MDL(rCAM), data.term = TRUE)

k <- 7

Aest <- Amat(rCAM, k)
Sest <- Smat(rCAM, k)

MGlist <- MGsforA(rCAM, 
                  K = k,
                  corner.strategy = 2) 

MG_tibble <- tibble(x = 1:k, y = MGlist[1:k]) %>%
  unnest(y) %>%
  inner_join(anno_mmus, by = c('y' = 'ensembl_gene_id'))

write_delim(MG_tibble, path = file.path(output_dir, "debCAM_MGs.txt"),  delim = "\t")

MGstat <- MGstatistic(data = data, #Calculates one versus everything fold change from subpop-specific expression profiles
                      A = Aest, 
                      boot.alpha = 0.05, 
                      nboot = 1000,
                      cores = 8)

MGlist.FC <- lapply(seq_len(k), function(x) 
  rownames(MGstat)[MGstat$idx == x & MGstat$OVE.FC > 10])

MGlist.FCboot <- lapply(seq_len(k), function(x) 
  rownames(MGstat)[MGstat$idx == x & MGstat$OVE.FC.alpha > 10])

MGlist.re <- reselectMG(data = data, 
                        MGlist = MGlist, 
                        fc.thres='q0.5', #default = q0.5
                        err.thres='q0.5') #default = NULL, q1.2 = the maximum reconstruction error times 1.2

rre <- redoASest(data = data, 
                 MGlist = MGlist.re, 
                 generalNMF = TRUE,
                 maxIter = 2, 
                 methy = FALSE)
#rre$Aest: re-estimated A matrix
#rre$Sest: re-estimated S matrix

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
simplexplot(data, Aest, MGlist, main=c('Initially detected markers'))
simplexplot(data, Aest, MGlist.FC, main=c('fc > 10'))
simplexplot(data, Aest, MGlist.FCboot, 
            main=c(expression(bold(paste('fc(bootstrap,',alpha,'=0.05) > 10')))))
simplexplot(data, Aest, MGlist.re, main=c("fc >= 'q0.2', error <= 'q1.2'"))

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
simplexplot(data, rre$Aest, MGlist, main=c('Initially detected markers'))
simplexplot(data, rre$Aest, MGlist.FC, main=c('fc > 10'))
simplexplot(data, rre$Aest, MGlist.FCboot, 
            main=c(expression(bold(paste('fc(bootstrap,',alpha,'=0.05) > 10')))))
simplexplot(data, rre$Aest, MGlist.re, main=c("fc >= 'q0.2', error <= 'q1.2'"))

