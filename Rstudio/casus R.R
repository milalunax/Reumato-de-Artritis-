#Deel1
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
align.Nor1 <- align(index = "ref_humaan", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2="SRR4785819_2_subset40k.fastq", output_file = "Nor1.BAM")
align.Nor2 <- align(index = "ref_humaan", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2="SRR4785820_2_subset40k.fastq", output_file = "Nor2.BAM")
align.Nor3 <- align(index = "ref_humaan", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2="SRR4785828_2_subset40k.fastq", output_file = "Nor3.BAM")
align.Nor4 <- align(index = "ref_humaan", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2="SRR4785831_2_subset40k.fastq", output_file = "Nor4.BAM")
align.RA1 <- align(index = "ref_humaan", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2="SRR4785979_2_subset40k.fastq", output_file = "RA1.BAM")
align.RA2 <- align(index = "ref_humaan", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2="SRR4785980_2_subset40k.fastq", output_file = "RA2.BAM")
align.RA3 <- align(index = "ref_humaan", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2="SRR4785986_2_subset40k.fastq", output_file = "RA3.BAM")
align.RA4 <- align(index = "ref_humaan", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2="SRR4785988_2_subset40k.fastq", output_file = "RA4.BAM")

library(Rsamtools)


#Deel 2
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
write.csv(countscasus, "RA_countmatrix.csv")

#Deel 3
casus_table=read.table("Reumato-de-Artritis-/count_matrix_RA.txt", row.names = 1, header = TRUE)
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
colnames(casus_table) <- rownames(treatment_casus_table)
dds_casus = DESeqDataSetFromMatrix(countData = casus_table,
                              colData = treatment_casus_table,
                              design = ~ treatment_casus)
dds_casus = DESeq(dds_casus)
resultaten <- results(dds_casus)
write.table(resultaten, file = 'Resultatencasus1.csv', row.names = TRUE, col.names = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
hoogste_fold_change
laagste_fold_change
laagste_p_waarde
#Eerst d volcano plot doen en daarna die GO analyse, daarna pas de pathway
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Reuma vs Normaal',
                pCutoff = 10e-32,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0, 
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                shape=c(1))
    
lab_italics <- paste0("italic('", rownames(resultaten), "')")
selectLab_italics = paste0(
  "italic('",
  c('VCAM1','KCTD12','ADAM12', 'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA'),
  "')")

EnhancedVolcano(resultaten,
                lab = lab_italics,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = selectLab_italics,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'pink', 'purple', 'red3'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
dev.copy(png, 'Volcanoplotcasus6.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()
#GO assay
install.packages("BiocManager")
BiocManager::install("goseq")
library(goseq)
BiocManager::install("geneLenDataBase")
library(geneLenDataBase)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(dplyr)
ALL= rownames(resultaten)
res <- as.data.frame(resultaten)
DEG= res%>%
  filter(padj<0.05)
DEG=rownames(DEG)


class(DEG)
gene.vector=as.integer(ALL%in%DEG)
names(gene.vector)=ALL 

head(gene.vector)
tail(gene.vector)


pwf=nullp(gene.vector,"hg19","geneSymbol")
GO.wall=goseq(pwf,"hg19","geneSymbol")


class(GO.wall)
head(GO.wall)
nrow(GO.wall)
enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)

library(dplyr)

top10 <- GO.wall %>%
  arrange(over_represented_pvalue) %>%   # kleinste p‑waarde = meest significant
  slice(1:10)
library(ggplot2)

ggplot(top10, aes(x = reorder(category, -over_represented_pvalue),
                  y = -log10(over_represented_pvalue))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 meest verrijkte GO‑categorieën",
    x = "GO‑categorie",
    y = "-log10(p‑waarde)"
  ) +
  theme_minimal()

ggplot(top10, aes(x = -log10(over_represented_pvalue),
                  y = reorder(category, over_represented_pvalue))) +
  geom_point(size = 4, color = "darkred") +
  labs(
    title = "Top 10 GO‑categorieën (dotplot)",
    x = "-log10(p‑waarde)",
    y = "GO‑categorie"
  ) +
  theme_minimal()

top10 <- top10 %>%
  mutate(ratio = numDEInCat / numInCat)

ggplot(top10, aes(
  x = ratio,
  y = reorder(category, ratio),
  size = numDEInCat,
  color = -log10(over_represented_pvalue)
)) +
  geom_point() +
  scale_color_viridis_c() +
  labs(
    title = "Top 10 GO categorieën (bubble plot)",
    x = "DE/total ratio",
    y = "GO categorie"
  ) +
  theme_minimal()

library(pheatmap)

mat <- -log10(as.matrix(top10$over_represented_pvalue))
rownames(mat) <- top10$category

pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
         color = viridis::viridis(50),
         main = "Top 10 GO categorieën")
res
res[1] <- NULL
res[2:5] <- NULL

#pathview(
  gene.data = res,
  pathway.id = "eco02026",  
  species = "eco",          
  gene.idtype = "KEGG",     
  limit = list(gene = 5)    
)