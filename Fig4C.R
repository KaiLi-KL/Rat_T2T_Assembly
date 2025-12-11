library(dplyr)
library(data.table)
library(vcfR)
library(rtracklayer)
library(GenomicRanges)
library(viridis)

chr_file         <- "genome.length.txt"
vcf_file         <- "rat_pg.svlen.vcf.gz"  
gene_gff_file    <- "gene.gff"                      
satellite_file  <- "SHRSP_Satellite.bed"


pdf_out          <- "per_chr_tracks_final.pdf"
per_chr_h        <- 0.4
left_label_pad   <- 1e7

sv_bin_size      <- 5e5
gene_bin_size    <- 5e5

rect_top_offset    <- 0.10
rect_bottom_offset <- -0.20
min_filled_len     <- 1e5
tick_by            <- 10e6
sv_palette <- colorRampPalette(viridis(256))

highlight_regions <- data.frame(
  CHROM = c("chr6","chr20","chr15"),
  start = c(138e6, 1e6, 29e6),
  end   = c(143e6, 6e6, 31.5e6),
  label = c("IGH","MHC","TCRa"),
  stringsAsFactors = FALSE
)

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
sv98 <- quantile(sv_bins$SV_count, 0.98, na.rm = TRUE)
sv99 <- quantile(sv_bins$SV_count, 0.99, na.rm = TRUE)
sv_max <- max(sv_bins$SV_count, na.rm = TRUE)

gene_gtf <- fread(gene_gff_file, sep = "\t", header = FALSE)
colnames(gene_gtf) <- c("CHROM","source","feature","start","end","score","strand","frame","attribute")
gene_bins <- gene_gtf %>%
  mutate(bin = floor(start / gene_bin_size)) %>%
  count(CHROM, bin, name = "gene_count") %>%
  mutate(bin_start = bin * gene_bin_size)
gene95 <- quantile(gene_bins$gene_count, 0.95, na.rm = TRUE)

satellite_bed <- fread(satellite_file, header = FALSE, col.names = c("CHROM","START","END"))
satellite <- satellite_bed  %>%
  mutate(CHROM = sub("^chr0*", "chr", CHROM))

chrom_layout <- data.frame(
  CHROM    = chromosomes,
  y_center = (n_chrom:1),
  stringsAsFactors = FALSE
)
chrom_layout$y0_box <- chrom_layout$y_center + rect_bottom_offset
chrom_layout$y1_box <- chrom_layout$y_center + rect_top_offset  

png("per_chr_tracks_final.png", 
    width = 4000, height = (n_chrom * per_chr_h + 1) * 300, res = 300)
par(mar = c(5, 6, 2, 4), xaxs = "i", yaxs = "i")

plot(0, type = "n",
     xlim = c(0, max(chrom_lengths$LENGTH)),
     ylim = c(0.5, n_chrom + 0.5),
     axes = FALSE, xlab = "", ylab = "")

n_col <- 100
sv_cols <- sv_palette(n_col)


for (i in seq_along(chromosomes)) {
  chr_name <- chromosomes[i]
  chr_y    <- chrom_layout$y_center[chrom_layout$CHROM == chr_name]
  chr_len  <- chrom_lengths$LENGTH[match(chr_name, chrom_lengths$CHROM)]
  
  y0_box <- chr_y + rect_bottom_offset
  y1_box <- chr_y + rect_top_offset
  
  axis_y <- y0_box - 0.05
  segments(0, axis_y, chr_len, axis_y)
  ticks <- seq(0, chr_len, by = tick_by)
  segments(ticks, axis_y - 0.02, ticks, axis_y + 0.02)
  labels <- ticks/1e6
  labels[1] <- "0"  
  text(ticks, rep(axis_y - 0.15, length(ticks)), 
       labels = labels, cex = 0.6, xpd = TRUE)

  track_height <- 0.2 
  
  chr_y0 <- chr_y + rect_bottom_offset +0.05 
  chr_y1 <- chr_y + rect_top_offset
  rect(0, chr_y0, chr_len, chr_y1, col = NA, border = "black")
  
  gene_y0 <- chr_y1 
  gene_y1 <- gene_y0 + track_height
  rect(0, gene_y0, chr_len, gene_y1, col = NA, border = "black")
 
  g_chr <- gene_bins %>% filter(CHROM == chr_name)
  if (nrow(g_chr)) {
    bar_top <- gene_y0 + pmin(g_chr$gene_count, gene95) / gene95 * (gene_y1 - gene_y0)
    rect(g_chr$bin_start, gene_y0,
         pmin(g_chr$bin_start + gene_bin_size, chr_len), bar_top,
         col = "#d95f02", border = NA)
  }
 
  s_chr <- satellite %>% filter(CHROM == chr_name)
  if (nrow(s_chr)) {
    gr <- GRanges(
      seqnames = Rle(s_chr$CHROM),
      ranges   = IRanges(start = s_chr$START, end = s_chr$END)
    )
    gr_merged <- reduce(gr, min.gapwidth = 1000)
    x0 <- start(gr_merged)
    x1 <- end(gr_merged)
    rect(x0, gene_y0, x1, gene_y1, col = "#ff99ff", border = NA)
  }

  gene_y0 <- chr_y1 
  gene_y1 <- gene_y0 + track_height
  rect(0, gene_y0, chr_len, gene_y1, col = NA, border = "black")
  
  sv_y0 <- chr_y + rect_bottom_offset +0.05
  sv_y1 <- chr_y + rect_top_offset
  rect(0, sv_y0, chr_len, sv_y1, col = NA, border = "black")

  s_chr <- sv_bins %>% filter(CHROM == chr_name)
  if (nrow(s_chr)) {
    val <- s_chr$SV_count
    idx <- pmax(1, ceiling(val / sv99 * n_col))
    rect(s_chr$bin_start, sv_y0,
         pmin(s_chr$bin_start + sv_bin_size, chr_len), sv_y1,
         col = sv_cols[idx], border = NA)
  }

  chr_y0 <- chr_y + rect_bottom_offset +0.05 
  chr_y1 <- chr_y + rect_top_offset
  rect(0, chr_y0, chr_len, chr_y1, col = NA, border = "black")

  text(-left_label_pad, chr_y, labels = chr_name, xpd = TRUE, adj = 1, font = 2)
  
}

if (nrow(highlight_regions)) {
  hi_rects <- merge(highlight_regions, chrom_layout, by = "CHROM", all.x = TRUE)
  hi_rects <- hi_rects[!is.na(hi_rects$y0_box), ]
  if (nrow(hi_rects)) {
    rect(hi_rects$start, hi_rects$y0_box,
         hi_rects$end,   hi_rects$y1_box,
         border = "red", col = NA, lwd = 2)
    text((hi_rects$start + hi_rects$end)/2 +1e6 ,
         hi_rects$y1_box + 0.35,
         labels = hi_rects$label,
         col = "blue", cex = 0.8, font = 1.5)
  }
}

legend(x = 2e8, y = 11.3,
       legend = c("Gene density", "Satellite", "Immune regions"),
       fill   = c("#d95f02", "#ff99ff", NA),
       border = c(NA, NA, "red"),
       bty = "n", cex = 0.8)

par(xpd = TRUE)
sv_breaks <- seq(0, sv99, length.out = n_col)
bar_x0 <- max(chrom_lengths$LENGTH) * 0.8
bar_x1 <- bar_x0 + 3e6
bar_y0 <- n_chrom - 15
bar_y1 <- n_chrom - 13

for (j in seq_along(sv_cols)) {
  rect(bar_x0,
       bar_y0 + (j-1)/n_col * (bar_y1 - bar_y0),
       bar_x1,
       bar_y0 + j/n_col * (bar_y1 - bar_y0),
       col = sv_cols[j], border = NA)
}

text(bar_x1 + 5e6, bar_y0 - 0.1, "0",  adj = 0, cex = 0.7)
text(bar_x1 + 5e6, bar_y1 + 0.1, sv99, adj = 0, cex = 0.7)
text((bar_x0+bar_x1)/2, bar_y1 + 0.5, "SV density", cex = 0.8)
dev.off()



##### Find co-localize region

# combined_bins <- sv_bins %>%
#   left_join(
#     gene_bins %>% select(CHROM, bin, gene_count),
#     by = c("CHROM", "bin")
#   ) %>%
#   mutate(gene_count = ifelse(is.na(gene_count), 0, gene_count))
# sv_threshold <- sv99
# gene_threshold <- gene95
# 
# high_density_regions <- combined_bins %>%
#   filter(
#     SV_count >= sv_threshold,
#     gene_count >= gene_threshold
#   ) %>%
#   mutate(bin_end = bin_start + sv_bin_size) %>%
#   select(CHROM, bin, bin_start, bin_end, SV_count, gene_count) %>%
#   arrange(CHROM, bin_start)
# 
# write.csv(x = high_density_regions , file = "high_density_regions.csv", row.names = FALSE,quote = FALSE)
# 
