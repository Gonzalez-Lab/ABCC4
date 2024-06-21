library(DESeq2)
library(dplyr)
library(pheatmap)
library(EnhancedVolcano)

### Culture BxPC3 mock vs MRP4+

counts.table <- read.csv2("PANC1 shMRP4 culture and BxPC3 MRP4+ culture xenograft.csv", header = TRUE) 

counts.table[which(counts.table$GeneSymbol=="CDH1"),]

count.matrixC <- dplyr::select(counts.table, c("mock.1","mock.2","mrp4.1","mrp4.2"))
colnames(count.matrixC) <- c('mock_1',"mock_2","MRP4_1","MRP4_2")
count.matrixC <- as.matrix(count.matrixC)

sample.table <- as.data.frame(matrix(c("mock","mock","MRP4","MRP4"), ncol = 1))
rownames(sample.table) <- c("mock_1","mock_2","MRP4_1","MRP4_2")
colnames(sample.table) <- c("group")

sample.table$group <- factor(sample.table$group)

levels(sample.table$group) 

##check tables are OK
all(rownames(sample.table) == colnames(count.matrixC))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
ddsC <- DESeqDataSetFromMatrix(countData=count.matrixC, 
                              colData=sample.table, 
                              design=~group)

# al objeto dds (que es el dataset) le agregamos los analisis
ddsC <- DESeq(ddsC)

resC <- results(ddsC)

# table(resC$padj != "NA")
table(resC$padj < 0.05)
table(resC$padj < 0.05 & resC$log2FoldChange< -1) # downregulados
table(resC$padj < 0.05 & resC$log2FoldChange> 1) # upregulados

#le agrego los genes a la tabla de resultados
resC$geneid <- counts.table$GeneSymbol

head(resC)

resC[which(resC$geneid=="ABCC4"),]

#A summary of the results can be generated:
summary(resC)

EnhancedVolcano(resC,
                lab = resC$geneid,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Culture BxPC3 mock vs MRP4',
                subtitle = "1055 DEGs (padj<0,05, log2FC>1)",
                selectLab = c('ABCC4'),
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 1.5,
                labSize = 5,
                col=c('black', 'gray', 'green', 'red3'),
                colAlpha = 1,
                xlim=c(-11,11))

#A PCA plot and a heatmap of the top genes: we need to use the rld object
rldC <- rlog(ddsC)

plotPCA(rldC, intgroup="group")

#A results table:
signif_padjC <- resC[which(resC$padj<0.05 & resC$padj!= "NA"),]
signif_padjC$UP_DOWN <- ifelse(signif_padjC$log2FoldChange>1,"UP",ifelse(signif_padjC$log2FoldChange< -1,"DOWN","..."))
head(signif_padjC)

signif.genes_C <- signif_padjC$geneid[which(signif_padjC$UP_DOWN!="...")] %>% na.omit()
signif.genes_CUP <- signif_padjC$geneid[which(signif_padjC$UP_DOWN=="UP")] %>% na.omit()
signif.genes_CDOWN <- signif_padjC$geneid[which(signif_padjC$UP_DOWN=="DOWN")] %>% na.omit()

#heatmap de DEGs
matC <- assay(rldC)[which(resC$geneid %in% signif.genes_C),]
matC <- matC - rowMeans(matC) # hago el z score al restar cada valor menos la media
head(matC)
pheatmap(matC, fontsize = 18)

write.csv(signif_padjC, file = "total signif genes culture BxPC3 mock vs MRP4+.csv")

write.csv2(signif_padjC, file = "total signif genes culture BxPC3 mock vs MRP4+.csv")

######################################################################################################################

library(enrichR)

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()

## if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", 
         "KEGG_2021_Human","MSigDB_Hallmark_2020","WikiPathway 2021 Human")

if (websiteLive) {
  enrichedC <- enrichr(c(signif.genes_CUP, signif.genes_CDOWN), dbs)
}

GO_BP.C <- enrichedC[["GO_Biological_Process_2023"]][which(enrichedC[["GO_Biological_Process_2023"]]$Adjusted.P.value<0.05),]
GO_BP.C$dataset <- rep("GO.BP", nrow(GO_BP.C))

GO_CC.C <- enrichedC[["GO_Cellular_Component_2023"]][which(enrichedC[["GO_Cellular_Component_2023"]]$Adjusted.P.value<0.05),]
GO_CC.C$dataset <- rep("GO.CC", nrow(GO_CC.C))

GO_MF.C <- enrichedC[["GO_Molecular_Function_2023"]][which(enrichedC[["GO_Molecular_Function_2023"]]$Adjusted.P.value<0.05),]
GO_MF.C$dataset <- rep("GO.MF", nrow(GO_MF.C))

MSigDB.C <- enrichedC[["MSigDB_Hallmark_2020"]][which(enrichedC[["MSigDB_Hallmark_2020"]]$Adjusted.P.value<0.05),]
MSigDB.C$dataset <- rep("MSigDB", nrow(MSigDB.C))

KEGG.C <- enrichedC[["KEGG_2021_Human"]][which(enrichedC[["KEGG_2021_Human"]]$Adjusted.P.value<0.05),]
KEGG.C$dataset <- rep("KEGG", nrow(KEGG.C))

Wiki.C <- enrichedC[["WikiPathway 2021 Human"]][which(enrichedC[["WikiPathway 2021 Human"]]$Adjusted.P.value<0.05),]
Wiki.C$dataset <- rep("Wiki", nrow(Wiki.C))

func.enrich.C <- rbind(GO_BP.C,GO_CC.C,GO_MF.C,MSigDB.C,KEGG.C,Wiki.C)

write.csv(func.enrich.C, "functional enrichment BxPC3 MRP4+ culture.csv")

write.csv2(func.enrich.C, "functional enrichment BxPC3 MRP4+ culture.csv")

#Plot Enrichr GO-BP output. (Plotting function contributed by I-Hsuan Lin)

if (websiteLive) {
  plotEnrich(enrichedC[[5]], showTerms = 5, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
}

