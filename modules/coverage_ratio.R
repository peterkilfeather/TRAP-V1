#Validation of Sudmant
coverage_ratio_sudmant <- read_delim(file.path(base_dir, "modules/coverage_ratios/coverage_ratios_sudmant.csv"), delim = ",") %>% 
  drop_na(ratio)

#Density plot facet
facet = T
facet_factor = "cell_type"

ggplot(coverage_ratio_sudmant, aes(ratio, fill = sample, colour = age)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_bw() +
  scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Density") +
  xlab("Termination codon ratio") +
  facet_grid(rows = vars(eval(as.name(facet_factor)))) + 
  theme(strip.text.y = element_text(angle = 360)) 

ggplot(coverage_ratio_sudmant, aes(ratio, fill = sample, colour = age)) + stat_ecdf(geom = "step") + 
  facet_grid(rows = vars(eval(as.name(facet_factor)))) + 
  theme_bw() + 
  scale_color_lancet() + 
  theme(strip.text.y = element_text(angle = 360)) +
  xlab("3\' UTR Termination Ratio") + 
  ylab("Cumulative Fraction")

ggsave(path = output_dir, filename="sudmant_cumulative_frac.png", width = width, height = height, dpi = dpi, units = units)

coverage_ratio_sudmant %>% 
  group_by(age, cell_type) %>% 
  summarise(mean_ratio = mean(ratio, na.rm = TRUE))

coverage_ratio_sudmant %>% 
  group_by(age, cell_type) %>% 
  filter(ratio >= 5) %>% 
  summarise(n())

d1_over_five <- coverage_ratio_sudmant %>% 
  group_by(age, cell_type) %>% 
  filter(ratio >= 5 & cell_type == "d1")

# genes <- unique(as.vector(d1_over_five$gene))
# 
# write.csv(genes, file = "over_five_ratio.csv")

#Application to PK1 and KW2
coverage_ratio_PK1_KW2 <- read_delim(file.path(base_dir, "modules/coverage_ratios/coverage_ratios_PK1_KW2.csv"), delim = ",") %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
  drop_na(ratio)

coverage_ratio_ip <- coverage_ratio_PK1_KW2 %>% 
  filter(ip == "ip")

coverage_ratio_ip_KW2 <- coverage_ratio_ip %>% 
  filter(region == "mb" & batch_trap == "b2" & ip_pool == "n")

#Density plot facet
facet = T
facet_factor_x = "region"
facet_factor_y = "ovx"

ggplot(coverage_ratio_ip, aes(ratio, fill = sample_name, colour = age)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_minimal_grid() +
  scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Density") +
  xlab("Termination codon ratio") +
  facet_grid(rows = vars(eval(as.name(facet_factor_y))))

ggplot(coverage_ratio_ip, aes(ratio, fill = sample_name, colour = age)) + 
  stat_ecdf(geom = "step") + 
  facet_grid(rows = vars(eval(as.name(facet_factor_x))), cols = vars(eval(as.name(facet_factor_y)))) + 
  theme_bw() + 
  scale_color_lancet() + 
  theme(strip.text.y = element_text(angle = 360)) +
  xlab("3\' UTR Termination Ratio") + 
  ylab("Cumulative Fraction")

ggsave(path = output_dir, filename = "pk_trap_termination_cumulative_frac.png", width = width, height = height, dpi = dpi, units = units)

ggplot(coverage_ratio_ip_KW2, aes(ratio, fill = sample_name, colour = age)) + 
  stat_ecdf(geom = "step") + 
  facet_grid(rows = vars(eval(as.name(facet_factor_y)))) + 
  theme_bw() + 
  scale_color_lancet()

mb_over_five_test <- coverage_ratio_ip_KW2 %>% 
  group_by(gene, age, ovx) %>% 
  summarise(mean_ratio = mean(ratio, na.rm = TRUE)) %>% 
  spread(ovx, mean_ratio) %>% 
  mutate(diff = ovx - wt)

over_five <- mb_over_five_test %>% 
  filter(abs(diff) > 1)

coverage_ratio_PK1_KW2 %>% 
  group_by(age, region, ip, ovx) %>% 
  summarise(n())

mb_over_five_old_ovx <- coverage_ratio_ip_KW2 %>% 
  filter(ratio >= 5 & region == "mb" & ip == "ip" & age == "old")
mb_over_five_old_ovx_genes <- unique(as.vector(mb_over_five_old_ovx$gene))
mb_over_five_young_ovx <- coverage_ratio_ip_KW2 %>% 
  filter(ratio >= 5 & region == "mb" & ip == "ip" & age == "young")
mb_over_five_young_ovx_genes <- unique(as.vector(mb_over_five_young_ovx$gene))

intersect(mb_over_five_old_ovx_genes, mb_over_five_young_ovx_genes)

old_not_young_genes <- setdiff(mb_over_five_old_ovx_genes, mb_over_five_young_ovx_genes)

young_not_old_genes <- setdiff(mb_over_five_young_ovx_genes, mb_over_five_old_ovx_genes)

coverage_ratio_ip_KW2 %>% 
  filter(gene %in% old_not_young_genes)
