# msa_sorted_colored_by_segments_sorted_brighter.R
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
library(ggplot2)
library(RColorBrewer)

msa_file <- "for_tree_cut.fas"
output_png <- "msa_sorted_colored_by_segments_sorted_brighter.png"

segments <- data.frame(
  start = c(1,1271,1562,1625,1968,2028,2811,3482,3647,3708,4393,4457,4514,
            6123,6184,6234,6271,6343,6409,6473,7503,7566,7636,7709),
  end   = c(1271,1561,1624,1967,2027,2811,3491,3646,3707,4392,4456,4523,6121,
            6183,6236,6270,6342,6408,6472,7493,7565,7635,7708,8131),
  label = c("ND4","ND4L","ARG-tRNA","ND3","GLY-tRNA","COX3","ATP6","ATP8","LYS-tRNA",
            "COX2","ASP-tRNA","SER-tRNA","COI","TYR-tRNA","CYS-tRNA",
            "putative_L_origin","ASN-tRNA","ALA-tRNA","TRP-tRNA","ND2",
            "MET-tRNA","GLN-tRNA","LEU-tRNA","CR2"),
  stringsAsFactors = FALSE
)

read_fasta_alignment <- function(path) {
  lines <- readLines(path, warn = FALSE)
  names <- character(); seqs <- character(); cur_name <- NULL; cur_seq <- NULL
  for (ln in lines) {
    if (nchar(ln) == 0) next
    if (substr(ln,1,1) == ">") {
      if (!is.null(cur_name)) { names <- c(names, cur_name); seqs <- c(seqs, paste0(cur_seq, collapse = "")) }
      cur_name <- sub("^>\\s*", "", ln); cur_seq <- character()
    } else { cur_seq <- c(cur_seq, gsub("\\s+", "", ln)) }
  }
  if (!is.null(cur_name)) { names <- c(names, cur_name); seqs <- c(seqs, paste0(cur_seq, collapse = "")) }
  names(seqs) <- names
  seqs
}

cat("Reading alignment:", msa_file, "\n")
seqs <- read_fasta_alignment(msa_file)
if (length(seqs) == 0) stop("No sequences found.")
lens <- unique(nchar(seqs))
if (length(lens) != 1) stop("Sequences have differing lengths in alignment.")
aln_len <- lens[1]
cat("Alignment length:", aln_len, "positions; sequences:", length(seqs), "\n")

extract_nums <- function(x) {
  m1 <- regmatches(x, regexpr("[0-9]+", x))
  main_num <- ifelse(lengths(m1) == 0, NA, as.numeric(m1))
  m2 <- regmatches(x, regexpr("_[0-9]+", x))
  sub_num <- ifelse(lengths(m2) == 0, NA, as.numeric(sub("_", "", m2)))
  data.frame(name = x, main_num = main_num, sub_num = sub_num, stringsAsFactors = FALSE)
}

mito_mask <- grepl("mitochondrion", names(seqs), ignore.case = TRUE)
non_mito_names <- names(seqs)[!mito_mask]
if (length(non_mito_names) > 0) {
  df_nums <- extract_nums(non_mito_names)
  df_nums$main_num_for_sort <- ifelse(is.na(df_nums$main_num), .Machine$integer.max - 1, df_nums$main_num)
  df_nums$sub_num_for_sort  <- ifelse(is.na(df_nums$sub_num), -1, df_nums$sub_num)
  df_nums <- df_nums[order(df_nums$main_num_for_sort, df_nums$sub_num_for_sort, df_nums$name), ]
  non_mito_sorted <- df_nums$name
} else { non_mito_sorted <- character(0) }

order_names <- c(non_mito_sorted, names(seqs)[mito_mask])
order_names <- order_names[order_names %in% names(seqs)]
seqs <- seqs[order_names]

cat("Sequences will be plotted in this order (top -> bottom):\n")
print(order_names)

is_gap <- function(ch) ch %in% c("-", ".")
aln_list <- vector("list", length(seqs))
i <- 1
for (nm in names(seqs)) {
  s <- strsplit(seqs[[nm]], "")[[1]]
  aln_list[[i]] <- data.frame(seq = nm,
                              pos = seq_along(s),
                              nuc = s,
                              is_gap = s %in% c("-", "."),
                              stringsAsFactors = FALSE)
  i <- i + 1
}
aln_df <- do.call(rbind, aln_list)

pos2gene <- rep(NA_character_, aln_len)
for (r in seq_len(nrow(segments))) {
  s <- segments$start[r]; e <- segments$end[r]; lab <- segments$label[r]
  if (s > aln_len) next
  e2 <- min(e, aln_len)
  pos2gene[s:e2] <- lab
}
aln_df$gene <- pos2gene[aln_df$pos]
aln_df$gene_fill <- ifelse(aln_df$is_gap, NA, aln_df$gene)

genes <- unique(segments$label)
n_genes <- length(genes)
if (n_genes <= 12) {
  base_pal <- brewer.pal(max(3, n_genes), "Paired") # яскравіші кольори
  pal <- colorRampPalette(base_pal)(n_genes)
} else {
  pal <- colorRampPalette(brewer.pal(12, "Paired"))(n_genes)
}
names(pal) <- genes

aln_df$seq <- factor(aln_df$seq, levels = rev(names(seqs)))

p <- ggplot(aln_df, aes(x = pos, y = seq)) +
  geom_tile(aes(fill = gene_fill), color = NA, height = 0.8) +  # трохи відступів між рядками
  scale_fill_manual(values = pal, na.value = "white", name = "Gene") +
  geom_vline(xintercept = seq(0, aln_len, by = 500), color = "gray90", size = 0.2) +
  geom_text(data = data.frame(seq = levels(aln_df$seq)),
            aes(x = aln_len + round(aln_len*0.02), y = seq,
                label = ifelse(grepl("R_", as.character(seq), ignore.case = TRUE), "←", "→")),
            inherit.aes = FALSE, size = 3) +
  labs(x = "Position", y = NULL,
       title = "Композиція NUMT",
       subtitle = "White = gap or position not in segments") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "right")

ggsave(output_png, p, width = 14, height = max(6, length(seqs) * 0.2), dpi = 300, units = "in")
cat("Saved:", output_png, "\n")

