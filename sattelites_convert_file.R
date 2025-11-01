##satellites for IGV
library(stringr)
library(Biostrings)

# читаємо fasta
fasta <- readDNAStringSet("sattelites_on_3_chr.fasta")

# заголовки
headers <- names(fasta)

# парсимо
seqid <- str_extract(headers,
                     "^[^:]+")
coords <- str_extract(headers,
                      "(?<=:)\\d+-\\d+")
start <- as.numeric(str_extract(coords,
                                "^[0-9]+"))
end   <- as.numeric(str_extract(coords, 
                                "(?<=-)\\d+"))
label <- str_trim(str_replace(headers, 
                              "^[^ ]+ ",
                              ""))

# формуємо GFF3-рядки
gff <- data.frame(
  seqid,
  source = "Manual",
  type = "satellite DRU-Sat-1",
  start,
  end,
  score = ".",
  strand = ".",
  phase = ".",
  attributes = paste0("ID=Seq",
                      seq_along(headers),
                      ";Name=", label),
  stringsAsFactors = FALSE
)

# пишемо у файл
cat(
  'track name="satellite DRU-Sat-1" description="Daboia russelii-like satellite DRU-Sat-1 sequence" visibility=2 type=gff3 useNameForDisplay=on\n',
  "##gff-version 3\n",
  "##sequence-region OX365966.1 1 212821320\n",
  file = "Vipera_sattelites.gff3"
)

write.table(
  gff, "Vipera_sattelites.gff3",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE
)







