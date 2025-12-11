library(dplyr)
library(data.table)
library(vcfR)
library(rtracklayer)
library(GenomicRanges)
library(viridis)

chr_file         <- "genome.length.txt"
vcf_file         <- "rat_pg.svlen.vcf.gz"  
gene_gff_file    <- "gene.gff" 

sv_bin_size   <- 1e5  
gene_bin_size <- 1e5 


ori_chrom_lengths <- fread(chr_file, header = FALSE,
                           col.names = c("CHROM","LENGTH"))
chrom_lengths <- ori_chrom_lengths %>%
  mutate(CHROM = sub("^chr0*", "chr", CHROM))
chromosomes <- chrom_lengths$CHROM
n_chrom <- length(chromosomes)


vcf <- read.vcfR(vcf_file)
fix_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
fix_df$POS   <- as.numeric(fix_df$POS)
gt <- vcf@gt[, -1]
gt_mat <- as.matrix(gt)
fix_df$count_sv <- apply(gt_mat, 1, function(x) sum(!is.na(x) & x != 0))



all_bins <- fix_df %>%
  group_by(CHROM) %>%
  summarise(
    max_bin = max(floor(POS / sv_bin_size), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  do(data.frame(CHROM = .$CHROM,
                bin = 0:.$max_bin)) %>%
  ungroup()


sv_counts_df <- fix_df %>%
  filter(count_sv > 0) %>%
  mutate(bin = floor(POS / sv_bin_size)) %>%
  group_by(CHROM, bin) %>%
  summarise(SV_count = sum(count_sv), .groups = "drop")


sv_bins <- all_bins %>%
  left_join(sv_counts_df, by = c("CHROM", "bin")) %>%
  mutate(
    SV_count = ifelse(is.na(SV_count), 0, SV_count),
    bin_start = bin * sv_bin_size
  )


sv95 <- quantile(sv_bins$SV_count, 0.95, na.rm = TRUE)
sv99 <- quantile(sv_bins$SV_count, 0.99, na.rm = TRUE)


gene_gtf <- fread(gene_gff_file, sep = "\t", header = FALSE)
colnames(gene_gtf) <- c("CHROM","source","feature","start","end","score","strand","frame","attribute")


gene_bins <- gene_gtf %>%
  mutate(bin = floor(start / gene_bin_size)) %>%
  count(CHROM, bin, name = "gene_count") %>%
  mutate(bin_start = bin * gene_bin_size)

gene95 <- quantile(gene_bins$gene_count, 0.95, na.rm = TRUE)


combined_bins <- all_bins %>%
  left_join(
    sv_counts_df %>% select(CHROM, bin, SV_count),
    by = c("CHROM", "bin")
  ) %>%
  left_join(
    gene_bins %>% select(CHROM, bin, gene_count),
    by = c("CHROM", "bin")
  ) %>%
  mutate(
    SV_count = ifelse(is.na(SV_count), 0, SV_count),
    gene_count = ifelse(is.na(gene_count), 0, gene_count)
  )


test_data <- combined_bins %>%
  filter(CHROM %in% chromosomes)

# Spearman
correlation_test_100kb <- cor.test(
  test_data$SV_count,
  test_data$gene_count,
  method = "spearman"
)


rho_value_100kb <- round(correlation_test_100kb$estimate, 3)
p_value_100kb <- format.pval(correlation_test_100kb$p.value, digits = 3, eps = 0.001)

# >= 95th Percentile
hotspot_data_100kb <- combined_bins %>%
  filter(SV_count >= sv95) %>%
  filter(CHROM %in% chromosomes)

hotspot_correlation_test <- cor.test(
  hotspot_data_100kb$SV_count,
  hotspot_data_100kb$gene_count,
  method = "spearman"
)

hotspot_rho_100kb <- round(hotspot_correlation_test$estimate, 3)
hotspot_p_100kb <- format.pval(hotspot_correlation_test$p.value, digits = 3, eps = 0.001)
