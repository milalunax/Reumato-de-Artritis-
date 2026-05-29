#Indexeren
setwd("C:/Users/mtous/OneDrive - NHL Stenden/transcriptomics/casus")
getwd()
install.packages('BiocManager')
BiocManager::install('Rsubread')
library(Rsubread)
buildindex(
  basename = 'ref_humaan',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)

#Mapping
align.Nor1 <- align(index = "ref_humaan", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2="SRR4785819_2_subset40k.fastq", output_file = "Nor1.BAM")
align.Nor2 <- align(index = "ref_humaan", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2="SRR4785820_2_subset40k.fastq", output_file = "Nor2.BAM")
align.Nor3 <- align(index = "ref_humaan", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2="SRR4785828_2_subset40k.fastq", output_file = "Nor3.BAM")
align.Nor4 <- align(index = "ref_humaan", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2="SRR4785831_2_subset40k.fastq", output_file = "Nor4.BAM")
align.RA1 <- align(index = "ref_humaan", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2="SRR4785979_2_subset40k.fastq", output_file = "RA1.BAM")
align.RA2 <- align(index = "ref_humaan", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2="SRR4785980_2_subset40k.fastq", output_file = "RA2.BAM")
align.RA3 <- align(index = "ref_humaan", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2="SRR4785986_2_subset40k.fastq", output_file = "RA3.BAM")
align.RA4 <- align(index = "ref_humaan", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2="SRR4785988_2_subset40k.fastq", output_file = "RA4.BAM")


#Countmatrix
count_matrix=featureCounts(
  files = "Nor1.BAM",
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE, 
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)
allcasussamples=c("Nor1.BAM", "Nor2.BAM", "Nor3.BAM", "Nor4.BAM", "RA1.BAM", "RA2.BAM", "RA3.BAM", "RA4.BAM")
count_matrix=featureCounts(
  files = allcasussamples,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE, 
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)
str(count_matrix)
countscasus <- count_matrix$counts
head(countscasus)
colnames(countscasus) <- c("Normaal1", "Normaal2.BAM", "Normaal3.BAM", "Normaal4.BAM", "Reuma1.BAM", "Reuma2.BAM", "Reuma3.BAM", "Reuma4.BAM")
head(countscasus)
write.csv(countscasus, "RheumatoidArthritis_countmatrix.csv")
View(read.csv("RheumatoidArthritis_subsetcountmatrix.csv"))

#Statistiek
  casus_table=read.table("count_matrix_RA.txt", row.names = 1, header=TRUE )
head(casus_table)

BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pathview")

library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)
treatment_casus <- c("Normaal", "Normaal", "Normaal", "Normaal", "Reuma", "Reuma", "Reuma", "Reuma")
treatment_casus_table <- data.frame(treatment_casus)
rownames(treatment_casus_table) <- c("Normaal1", "Normaal2", "Normaal3", "Normaal4", "Reuma1", "Reuma2", "Reuma3", "Reuma4")
head(treatment_casus_table)
colnames(casus_table)
rownames(treatment_casus_table)
ncol(casus_table)
colnames(casus_table) <- rownames(treatment_casus_table)
dds_casus = DESeqDataSetFromMatrix(countData = casus_table,
                                   colData = treatment_casus_table,
                                   design = ~ treatment_casus)
dds_casus = DESeq(dds_casus)
resultaten <- results(dds_casus)
write.table(resultaten, file = 'DESeq2_resultaten.csv', row.names = TRUE, col.names = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
hoogste_fold_change
laagste_fold_change
laagste_p_waarde
resultaten

#Volcano plot
library(EnhancedVolcano)
topgenes <- head(rownames(resultaten[order(resultaten$padj, na.last = NA), ]), 10)

EnhancedVolcano(
  resultaten,
  lab = rownames(resultaten),
  x = 'log2FoldChange',
  y = 'pvalue',
  
  title = 'De differentiële genexpressie (DESeq2) ',
  subtitle = 'RA-samples vergeleken controle-samples',
  
  col = c(
    'grey80',      # NS
    '#A8DADC',     # pastel blauwgroen (alleen p)
    '#F6BD60',     # pastel geel/oranje (alleen FC)
    '#F28482'      # pastel zalm/rood (p + FC)
  ),
  
  pCutoff = 0.05,
  FCcutoff = 1,
  
  selectLab = topgenes,
  
  pointSize = 2.5,
  labSize = 4.5,
  labCol = '#D62828',
  labFace = 'bold',
  
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 12,
  
  drawConnectors = TRUE,
  widthConnectors = 0.8,
  colConnectors = '#2A9D8F',
  
  xlab = "Log2 fold change",
  ylab = "-Log10 p-value"
)

dev.copy(png, 'RheumatoidArthritisVolcanoplot.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

#GO-analyse
install.packages("BiocManager")
BiocManager::install("goseq")
library(goseq)
BiocManager::install("geneLenDataBase")
library(geneLenDataBase)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(dplyr)
library("magrittr")

ALL= rownames(resultaten)
res <- as.data.frame(resultaten)
DEG= res%>%
  filter(padj<0.05, log2FoldChange<+-1 | log2FoldChange>=1)

DEG=rownames(DEG)
DEG
res
class(DEG)
gene.vector=as.integer(ALL%in%DEG)
names(gene.vector)=ALL 

head(gene.vector)
tail(gene.vector)

pwf=nullp(gene.vector,"hg19","geneSymbol")
plot(pwf)
GO.wall=goseq(pwf,"hg19","geneSymbol")

dev.copy(png, 'RheumatoidArthritisPWF.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

class(GO.wall)
head(GO.wall)
nrow(GO.wall)
enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)

library(ggplot2)

top10 <- GO.wall %>%
  arrange(over_represented_pvalue) %>%   # kleinste p‑waarde = meest significant
  slice(1:10)

GO.wall$diff <- GO.wall$numDEInCat / GO.wall$numInCat
GO.wall$abs_diff <- GO.wall$numInCat - GO.wall$numDEInCat
GO.wall$ratio <- GO.wall$numDEInCat / GO.wall$numInCat

top10 <- GO.wall %>%
  arrange(over_represented_pvalue) %>%
  slice(1:10)
top10$negLogP=-log10(top10$over_represented_pvalue)

top10$ratio_pct <- (top10$numDEInCat / top10$numInCat) * 100

top10 <- top10 %>%
  mutate(ratio_pct = ratio * 100)
ggplot(top10, aes(
  x = ratio_pct,
  y = reorder(term, ratio_pct)
)) +
  geom_col(fill = "#2C7BB6") +
  labs(
    title = "Top 10 genen van oververtegenwoordigde biologische processen",
    x = "Percentage DE-genen",
    y = "Biologisch proces"
  ) +
  theme_minimal(base_size = 14)
ggsave("RheumatoidArthritis-GO-analyse.png", width = 15, height = 8, dpi = 300)

#pathway-analyse
res
res[1] <- NULL
res[2:5] <- NULL

pathview(
  gene.data = res,
  pathway.id = "hsa04062",  
  species = "hsa",          
  gene.idtype = "SYMBOL",     
  limit = list(gene = 5)    
)
