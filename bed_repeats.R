## Конвертація gff3 у bed для IGV
## install.packages("dplyr")   
## якщо ще не встановлено
library(dplyr)
library(stringr)
getwd()

## Файл
gff_file <- "Galaxy55-[RepeatMasker repeat annotation on data 8].gff"

## Зчитування GFF
gff <- read.delim(gff_file, 
                  header = FALSE,
                  comment.char = "#",
                  stringsAsFactors = FALSE)

## Колонки GFF: 1=chrom, 2=source, 3=type, 4=start, 
##5=end, 6=score, 7=strand, 8=phase, 9=attributes
colnames(gff) <- c("chrom",
                   "source",
                   "type",
                   "start",
                   "end",
                   "score",
                   "strand",
                   "phase",
                   "attributes")

## Motif та Class з attributes
gff <- gff %>%
  rowwise() %>%
  mutate(
    motif = str_extract(attributes, 'Motif:[^;]+'),
    class = str_extract(attributes, 'Class=[^;]+'),
    motif = str_replace_all(motif, 'Motif:|"', ''),   
    class = str_replace_all(class, 'Class=|"', ''),   
    name = paste0(motif, "_", class)                  ## комбінований підпис
  )

##  BED-подібний датафрейм
bed <- gff %>%
  select(chrom, start, end, name) %>%
  mutate(start = start - 1)   # 0-based start для BED

# файл
write.table(bed, "repeats.bed", 
            sep="\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
