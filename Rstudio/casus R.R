setwd("C:/Users/mtous/OneDrive - NHL Stenden/transcriptomics/casus/")
getwd()
install.packages('BiocManager')
BiocManager::install('Rsubread')
library(Rsubread)
buildindex(
  basename = 'ref_humaan',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)
align.fw.RA1 <- align(index = "ref_humaan", readfile1 = "", output_file = "eth1.BAM")