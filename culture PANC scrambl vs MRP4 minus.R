library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)

### Culture PANC1 scrambl vs shot hairpin MRP4

counts.table <- read.csv2("PANC1 shMRP4 culture and BxPC3 MRP4+ culture xenograft.csv", header = TRUE) 

counts.table[which(counts.table$GeneSymbol=="CDH1"),]

count.matrixP <- dplyr::select(counts.table, c("SCR.1","SCR.2","SH.1","SH.2"))
colnames(count.matrixP) <- c('scrambl_1',"scrambl_2","MRP4sh_1","MRP4sh_2")
count.matrixP <- as.matrix(count.matrixP)

sample.table <- as.data.frame(matrix(c("scrambl","scrambl","MRP4sh","MRP4sh"), ncol = 1))
rownames(sample.table) <- c("scrambl_1","scrambl_2","MRP4sh_1","MRP4sh_2")
colnames(sample.table) <- c("group")

sample.table$group <- factor(sample.table$group)

levels(sample.table$group)
sample.table$group <- relevel(sample.table$group, "scrambl")

#check tables are OK
all(rownames(sample.table) == colnames(count.matrixP))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
ddsP <- DESeqDataSetFromMatrix(countData=count.matrixP, 
                               colData=sample.table, 
                               design=~group)

# al objeto dds (que es el dataset) le agregamos los analisis
ddsP <- DESeq(ddsP)

resP <- results(ddsP)

table(resP$padj != "NA")
table(resP$padj < 0.05)
table(resP$padj < 0.05 & resP$log2FoldChange  < -1) # downregulated
table(resP$padj < 0.05 & resP$log2FoldChange > 1) # upregulated

#add gene symbol to results table
resP$geneid <- counts.table$GeneSymbol

head(resP)

gene <- "GATA2"
resP[which(resP$geneid==gene),]

#A summary of the results can be generated:
summary(resP)

EnhancedVolcano(resP,
                lab = resP$geneid,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'PANC scrambl vs shMRP4',
                subtitle = "212 DEGs (padj<0,05, log2FC>1)",
                selectLab = c('ABCC4'),
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 1.5,
                labSize = 5,
                col=c('black', 'gray', 'green', 'red3'),
                colAlpha = 1,
                xlim=c(-10,10))

#A PCA plot and a heatmap of the top genes: we need to use the rld object
rldP <- rlog(ddsP)

plotPCA(rldP, intgroup="group")

#A results table:
signif_padjP <- resP[which(resP$padj<0.05 & resP$padj!= "NA"),]
signif_padjP$UP_DOWN <- ifelse(signif_padjP$log2FoldChange>1,"UP",ifelse(signif_padjP$log2FoldChange< -1,"DOWN","..."))
head(signif_padjP)

signif.genes_P <- signif_padjP$geneid[which(signif_padjP$UP_DOWN!="...")] %>% na.omit()
signif.genes_PUP <- signif_padjP$geneid[which(signif_padjP$UP_DOWN=="UP")] %>% na.omit()
signif.genes_PDOWN <- signif_padjP$geneid[which(signif_padjP$UP_DOWN=="DOWN")] %>% na.omit()

#heatmap de DEGs
matP <- assay(rldP)[which(resP$geneid %in% signif.genes_P),]
matP <- matP - rowMeans(matP) # hago el z score al restar cada valor menos la media
head(matP)
pheatmap(matP, fontsize = 18)

write.csv2(signif_padjP, file = "total signif genes PANC scrambl vs shMRP4.csv")

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
  enrichedP <- enrichr(c(signif.genes_PUP, signif.genes_PDOWN), dbs)
}


GO_BP.P <- enrichedP[["GO_Biological_Process_2023"]][which(enrichedP[["GO_Biological_Process_2023"]]$Adjusted.P.value<0.05),]
GO_BP.P$dataset <- rep("GO.BP", nrow(GO_BP.P))

GO_CC.P <- enrichedP[["GO_Cellular_Component_2023"]][which(enrichedP[["GO_Cellular_Component_2023"]]$Adjusted.P.value<0.05),]
GO_CC.P$dataset <- rep("GO.CC", nrow(GO_CC.P))

GO_MF.P <- enrichedP[["GO_Molecular_Function_2023"]][which(enrichedP[["GO_Molecular_Function_2023"]]$Adjusted.P.value<0.05),]
GO_MF.P$dataset <- rep("GO.MF", nrow(GO_MF.P))

MSigDB.P <- enrichedP[["MSigDB_Hallmark_2020"]][which(enrichedP[["MSigDB_Hallmark_2020"]]$Adjusted.P.value<0.05),]
MSigDB.P$dataset <- rep("MSigDB", nrow(MSigDB.P))

KEGG.P <- enrichedP[["KEGG_2021_Human"]][which(enrichedP[["KEGG_2021_Human"]]$Adjusted.P.value<0.05),]
KEGG.P$dataset <- rep("KEGG", nrow(KEGG.P))

Wiki.P <- enrichedP[["WikiPathway 2021 Human"]][which(enrichedP[["WikiPathway 2021 Human"]]$Adjusted.P.value<0.05),]
Wiki.P$dataset <- rep("Wiki", nrow(Wiki.P))

func.enrich.P <- rbind(GO_BP.P,GO_CC.P,GO_MF.P,MSigDB.P,KEGG.P,Wiki.P)

write.csv(func.enrich.P, "functional enrichment PANC1 MRP4-.csv")

write.csv2(func.enrich.P, "functional enrichment PANC1 MRP4-.csv")
  
#Plot Enrichr GO-BP output. (Plotting function contributed by I-Hsuan Lin)
  
if (websiteLive) {
    plotEnrich(enrichedP[[5]], showTerms = 5, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
  }

