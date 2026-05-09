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
align.Nor1 <- align(index = "ref_humaan", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2="SRR4785819_2_subset40k.fastq", output_file = "Nor1.BAM")
align.Nor2 <- align(index = "ref_humaan", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2="SRR4785820_2_subset40k.fastq", output_file = "Nor2.BAM")
align.Nor3 <- align(index = "ref_humaan", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2="SRR4785828_2_subset40k.fastq", output_file = "Nor3.BAM")
align.Nor4 <- align(index = "ref_humaan", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2="SRR4785831_2_subset40k.fastq", output_file = "Nor4.BAM")
align.RA1 <- align(index = "ref_humaan", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2="SRR4785979_2_subset40k.fastq", output_file = "RA1.BAM")
align.RA2 <- align(index = "ref_humaan", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2="SRR4785980_2_subset40k.fastq", output_file = "RA2.BAM")
align.RA3 <- align(index = "ref_humaan", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2="SRR4785986_2_subset40k.fastq", output_file = "RA3.BAM")
align.RA4 <- align(index = "ref_humaan", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2="SRR4785988_2_subset40k.fastq", output_file = "RA4.BAM")
