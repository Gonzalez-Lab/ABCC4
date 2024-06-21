library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)

### Xenograft BxPC3 mock vs MRP4+

counts.table <- read.csv2("PANC1 shMRP4 culture and BxPC3 MRP4+ culture xenograft.csv", header = TRUE) 

counts.table[which(counts.table$GeneSymbol=="CDH1"),]

count.matrixT <- dplyr::select(counts.table, c("Xmock.1","Xmock.2","Xmock.3","Xmrp4.1","Xmrp4.2","Xmrp4.3"))
colnames(count.matrixT) <- c("mock_1","mock_2","mock_3","MRP4_1", "MRP4_2","MRP4_3")
count.matrixT <- as.matrix(count.matrixT)

sample.tableT <- as.data.frame(matrix(c("mock","mock","mock","MRP4","MRP4","MRP4"), ncol = 1))
rownames(sample.tableT) <- c("mock_1","mock_2","mock_3","MRP4_1", "MRP4_2","MRP4_3")
colnames(sample.tableT) <- c("group")

sample.tableT$group <- factor(sample.tableT$group)

levels(sample.tableT$group) 

#chequeo que las tablas estan bien
all(rownames(sample.tableT) == colnames(count.matrixT))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
ddsT <- DESeqDataSetFromMatrix(countData=count.matrixT, 
                               colData=sample.tableT, 
                               design=~group)

# al objeto dds (que es el dataset) le agregamos los analisis
ddsT <- DESeq(ddsT)

resT <- results(ddsT)

table(resT$padj != "NA")
table(resT$padj < 0.05)
table(resT$padj < 0.05 & resT$log2FoldChange< -1) # downregulated
table(resT$padj < 0.05 & resT$log2FoldChange> 1) # upregulated

#add gene symbol to results table
resT$geneid <- counts.table$GeneSymbol

head(resT)

resT[which(resT$geneid=="ABCC4"),]

#A summary of the results can be generated:
summary(resT)

EnhancedVolcano(resT,
                lab = resT$geneid,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Xenograft BxPC3 mock vs MRP4',
                subtitle = "1254 DEG (FDR<0,05, log2FC>1)",
                selectLab = c('ABCC4'),
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 1.5,
                labSize = 5,
                col=c('black', 'gray', 'green', 'red3'),
                colAlpha = 1,
                xlim=c(-15,15))

#A PCA plot and a heatmap of the top genes: we need to use the rld object
rldT <- rlog(ddsT)

plotPCA(rldT, intgroup="group")


#A results table:
signif_padjT <- resT[which(resT$padj<0.05 & resT$padj!= "NA"),]
signif_padjT$UP_DOWN <- ifelse(signif_padjT$log2FoldChange>1,"UP",ifelse(signif_padjT$log2FoldChange< -1,"DOWN","..."))
head(signif_padjT)

signif.genes_T <- signif_padjT$geneid[which(signif_padjT$UP_DOWN!="...")] %>% na.omit()
signif.genes_TUP <- signif_padjT$geneid[which(signif_padjT$UP_DOWN=="UP")] %>% na.omit()
signif.genes_TDOWN <- signif_padjT$geneid[which(signif_padjT$UP_DOWN=="DOWN")] %>% na.omit()

#heatmap de DEGs
matT <- assay(rldT)[which(resT$geneid %in% signif.genes_T),]
matT <- matT - rowMeans(matT) # hago el z score al restar cada valor menos la media
head(matT)
pheatmap(matT, fontsize = 18)

write.csv2(signif_padjT, file = "total signif genes xenograft BxPC3 mock vs MRP4+.csv")

#############################################################################################################
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
  enrichedT <- enrichr(c(signif.genes_TUP, signif.genes_TDOWN), dbs)
}

GO_BP.T <- enrichedT[["GO_Biological_Process_2023"]][which(enrichedT[["GO_Biological_Process_2023"]]$Adjusted.P.value<0.05),]
GO_BP.T$dataset <- rep("GO.BP", nrow(GO_BP.T))

GO_CC.T <- enrichedT[["GO_Cellular_Component_2023"]][which(enrichedT[["GO_Cellular_Component_2023"]]$Adjusted.P.value<0.05),]
GO_CC.T$dataset <- rep("GO.CC", nrow(GO_CC.T))

GO_MF.T <- enrichedT[["GO_Molecular_Function_2023"]][which(enrichedT[["GO_Molecular_Function_2023"]]$Adjusted.P.value<0.05),]
GO_MF.T$dataset <- rep("GO.MF", nrow(GO_MF.T))

MSigDB.T <- enrichedT[["MSigDB_Hallmark_2020"]][which(enrichedT[["MSigDB_Hallmark_2020"]]$Adjusted.P.value<0.05),]
MSigDB.T$dataset <- rep("MSigDB", nrow(MSigDB.T))

KEGG.T <- enrichedT[["KEGG_2021_Human"]][which(enrichedT[["KEGG_2021_Human"]]$Adjusted.P.value<0.05),]
KEGG.T$dataset <- rep("KEGG", nrow(KEGG.T))

Wiki.T <- enrichedT[["WikiPathway 2021 Human"]][which(enrichedT[["WikiPathway 2021 Human"]]$Adjusted.P.value<0.05),]
Wiki.T$dataset <- rep("Wiki", nrow(Wiki.T))

func.enrich.T <- rbind(GO_BP.T,GO_CC.T,GO_MF.T,MSigDB.T,KEGG.T,Wiki.T)

write.csv(func.enrich.T, "functional enrichment BxPC3 MRP4+ tumor.csv")

write.csv2(func.enrich.T, "functional enrichment BxPC3 MRP4+ tumor.csv")

#Plot Enrichr GO-BP output. (Plotting function contributed by I-Hsuan Lin)

if (websiteLive) {
  plotEnrich(enrichedT[[5]], showTerms = 5, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
}

