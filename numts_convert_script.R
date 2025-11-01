##create a NUMT gff3
library(dplyr)

##segment
segments <- data.frame(
  start = c(102521842,
            102461295,
            102608636,
            102634743, 
            102695707,
            102408254,
            102468996,
            102529554,
            102616750,
            102642210,
            102669701,
            102403203, 
            102491632,
            102664339,
            102480857,
            102480225),
  end = c(102529972,
          102469414, 
          102611817,
          102637276, 
          102697530,
          102411518, 
          102469414, 
          102529972,
          102617706,
          102643166,
          102670641,
          102408252,
          102495635,
          102665001, 
          102484742, 
          102480492),
  label = c("8131",
            "8119",
            "3181",
            "2533", 
            "1823",
            "3264",
            "418",
            "418",
            "956",
            "956",
            "940", 
            "5049", 
            "4003", 
            "662", 
            "3885",
            "267")
)

## GFF3
gff_numt <- segments %>%
  rowwise() %>%
  mutate(
    chrom = "chr3",              #
    source = "Manual",
    type = "NUMT",
    score = ".",
    strand = ".",                
    phase = ".",
    attributes = paste0("ID=NUMT_",
                        label, 
                        ";Name=NUMT_",
                        label)
  ) %>%
  select(chrom,
         source, 
         type,
         start,
         end, 
         score,
         strand,
         phase, 
         attributes)

# Добавляем track line для IGV
track_line <- "track name=\"NUMTs\" description=\"NUMT pseudogenes\" visibility=2 type=gff3 useNameForDisplay=on"

file_conn <- file("numts_simple.gff3")
writeLines(track_line, 
           file_conn)
close(file_conn)

# Записываем данные
write.table(gff_numt, 
            "numts_simple.gff3",
            sep="\t", 
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE,
            append=TRUE)

