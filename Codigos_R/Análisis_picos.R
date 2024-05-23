## Análisis de ChIP con ChIPseeker
##Elaborado por Luis Rivera
## Inmunoprecipitación de cromatina con anti-p65 y anti-Nrf2 en línea SH-SY5Y

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(rtracklayer)
library(ReactomePA)
library(consensusSeekeR)
library(GenomeInfoDb)
library(ggVennDiagram)
library(ggplot2)
library(tibble)

#se cargan los archivos .bed de las réplicas unidas. En este caso serían las
#réplicas que según el IDR generan las mejores parejas. La unión se da usando
#bedtools intersect.

peak_CN_p65 <-readPeakFile(file.path("CN_p65_14-idr.narrowPeak"), as="GRanges")
peak_CP_p65 <-readPeakFile(file.path("CP_p65_24-idr.narrowPeak"), as="GRanges")
peak_CP_Nrf2 <-readPeakFile(file.path("CP_Nrf2-idr_23.narrowPeak"), as="GRanges")
peak_Resv_p65 <-readPeakFile(file.path("Resv_p65_34-idr.narrowPeak"), as="GRanges")
peak_Resv_Nrf2 <-readPeakFile(file.path("Resv_Nrf2-idr_34.narrowPeak"), as="GRanges")
peak_1422_p65 <-readPeakFile(file.path("1422_p65_24-idr.narrowPeak"), as="GRanges")
peak_1422_Nrf2 <-readPeakFile(file.path("1422_Nrf2-idr_24.narrowPeak"), as="GRanges")

files <- getSampleFiles()

file_p65_resv <- list(
  CN = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/CN_p65_14-idr.narrowPeak",
  CP = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/CP_p65_24-idr.narrowPeak",
  Resv = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/Resv_p65_34-idr.narrowPeak")

file_p65_SWitA <- list(
  CN = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/CN_p65_14-idr.narrowPeak",
  CP = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/CP_p65_24-idr.narrowPeak",
  SWitA = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/1422_p65_24-idr.narrowPeak")


file_Nrf2 <- list(
  CP = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/CP_Nrf2-idr_23.narrowPeak", 
  Resv = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/Resv_Nrf2-idr_34.narrowPeak", 
  SWitA = "C:/Users/Luis/OneDrive/4. CIBCM/Resultados ensayos/ChIP/IDR_picos/1422_Nrf2-idr_24.narrowPeak")


#ChIP peaks coverage plot

covplot(peak_CN_p65, weightCol="V7", xlab = "Tama�o del cromosoma (pb)" , ylab = "", title = "Condici�n:CN Anticuerpo:p65")
covplot(peak_CP_p65, weightCol="V7", xlab = "Tama�o del cromosoma (pb)" , ylab = "", title = "Condici�n:CP Anticuerpo:p65")
covplot(peak_CP_Nrf2, weightCol="V7", xlab = "Tama�o del cromosoma (pb)" , ylab = "", title = "Condici�n:CP Anticuerpo:Nrf2")
covplot(peak_Resv_p65, weightCol="V7", xlab = "Tama�o del cromosoma (pb)" , ylab = "", title = "Condici�n:Resveratrol 20 �M Anticuerpo:p65")
covplot(peak_Resv_Nrf2, weightCol="V7", xlab = "Tama�o del cromosoma (pb)" , ylab = "", title = "Condici�n:Resveratrol 20 �M Anticuerpo:Nrf2")
covplot(peak_1422_p65, weightCol="V7", xlab = "Tama�o del cromosoma (pb)" , ylab = "", title = "Condici�n:SWitA 0.5 �M Anticuerpo:p65")
covplot(peak_1422_Nrf2, weightCol="V7", xlab = "Tama�o del cromosoma (pb)" , ylab = "", title = "Condici�n:SWitA 0.5 �M Anticuerpo:Nrf2")


#Profile of ChIP peaks binding to TSS regions

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_CN_p65 <- getTagMatrix(peak_CN_p65, windows=promoter)
tagMatrix_CP_p65 <- getTagMatrix(peak_CP_p65, windows=promoter)
tagMatrix_CP_Nrf2 <- getTagMatrix(peak_CP_Nrf2, windows=promoter)
tagMatrix_Resv_p65 <- getTagMatrix(peak_Resv_p65, windows=promoter)
tagMatrix_Resv_Nrf2 <- getTagMatrix(peak_Resv_Nrf2, windows=promoter)
tagMatrix_SWitA_p65 <- getTagMatrix(peak_1422_p65, windows=promoter)
tagMatrix_SWitA_Nrf2 <- getTagMatrix(peak_1422_Nrf2, windows=promoter)


tagHeatmap(tagMatrix_CN_p65, title = "Condici�n:CN Anticuerpo:p65") #Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix_CP_p65,title = "Condici�n:CP Anticuerpo:p65") #Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix_CP_Nrf2,title = "Condici�n:CP Anticuerpo:Nrf2") #Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix_Resv_p65,title = "Condici�n:Resv 20 �M Anticuerpo:p65") #Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix_Resv_Nrf2,title = "Condici�n:Resv 20 �M Anticuerpo:Nrf2") #Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix_SWitA_p65,title = "Condici�n:SWitA 0.5 �M Anticuerpo:p65") #Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix_SWitA_Nrf2,title = "Condici�n:SWitA 0.5 �M Anticuerpo:Nrf2") #Heatmap of ChIP binding to TSS regions


#Heatmaps. Se debe cambiar cada dataset.

peakHeatmap(peak = peak_1422_p65,
            TxDb = txdb,
            upstream = rel(0.2),
            downstream = rel(0.2),
            by = "gene",
            type = "body",
            nbin = 800)


#Average Profile of ChIP peaks binding to TSS region
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix_CN_p65, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
plotAvgProf(tagMatrix_CP_p65, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
plotAvgProf(tagMatrix_CP_Nrf2, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
plotAvgProf(tagMatrix_Resv_p65, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
plotAvgProf(tagMatrix_Resv_Nrf2, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
plotAvgProf(tagMatrix_SWitA_p65, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
plotAvgProf(tagMatrix_SWitA_Nrf2, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

plotPeakProf2(peak = peak_1422_p65, upstream = 3000, downstream = 3000,
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V7",ignore_strand = F)

plotPeakProf2(peak = peak_1422_p65, upstream = 3000, downstream = 3000,
              conf = 0.95, by = "gene", type = "start_site", nbin = 800,
              TxDb = txdb, weightCol = "V7",ignore_strand = F)


TSS_matrix <- getTagMatrix(peak = peak_CN_p65, 
                              TxDb = txdb,
                              upstream = 3000,
                              downstream = 3000, 
                              type = "end_site",
                              by = "gene",
                              weightCol = "V5",
                              nbin = 50)

plotPeakProf(tagMatrix = TSS_matrix, conf = 0.95)

################################################################

#Peak Annotation

peakAnno_CN_p65 <- annotatePeak(peak_CN_p65, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_CP_p65 <- annotatePeak(peak_CP_p65, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_CP_Nrf2 <- annotatePeak(peak_CP_Nrf2, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_Resv_p65 <- annotatePeak(peak_Resv_p65, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_Resv_Nrf2 <- annotatePeak(peak_Resv_Nrf2, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_SWitA_p65 <- annotatePeak(peak_1422_p65, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_SWitA_Nrf2 <- annotatePeak(peak_1422_Nrf2, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")

plotAnnoPie(peakAnno_CN_p65)
plotAnnoPie(peakAnno_CP_p65)
plotAnnoPie(peakAnno_CP_Nrf2)
plotAnnoPie(peakAnno_Resv_p65)
plotAnnoPie(peakAnno_Resv_Nrf2)
plotAnnoPie(peakAnno_SWitA_p65)
plotAnnoPie(peakAnno_SWitA_Nrf2)


plotDistToTSS(peakAnno_SWitA_p65,
              title="Distribution of transcription factor-binding loci relative to TSS")



#################################################################
#Functional enrichment analysis

pathway_CN_p65 <- enrichPathway(as.data.frame(peakAnno_CN_p65)$geneId)
pathway_name_CN_p65<- setReadable(pathway_CN_p65, 'org.Hs.eg.db', 'ENTREZID') 

pathway_CP_p65 <- enrichPathway(as.data.frame(peakAnno_CP_p65)$geneId)
pathway_name_CP_p65<- setReadable(pathway_CP_p65, 'org.Hs.eg.db', 'ENTREZID') 

pathway_CP_Nrf2 <- enrichPathway(as.data.frame(peakAnno_CP_Nrf2)$geneId)
pathway_name_CP_Nrf2<- setReadable(pathway_CP_Nrf2, 'org.Hs.eg.db', 'ENTREZID') 

pathway_Resv_p65 <- enrichPathway(as.data.frame(peakAnno_Resv_p65)$geneId)
pathway_name_Resv_p65<- setReadable(pathway_Resv_p65, 'org.Hs.eg.db', 'ENTREZID') 

pathway_Resv_Nrf2 <- enrichPathway(as.data.frame(peakAnno_Resv_Nrf2)$geneId)
pathway_name_Resv_Nrf2<- setReadable(pathway_Resv_Nrf2, 'org.Hs.eg.db', 'ENTREZID') 

pathway_SWitA_p65 <- enrichPathway(as.data.frame(peakAnno_SWitA_p65)$geneId)
pathway_name_SWitA_p65<- setReadable(pathway_SWitA_p65, 'org.Hs.eg.db', 'ENTREZID') 

pathway_SWitA_Nrf2 <- enrichPathway(as.data.frame(peakAnno_SWitA_Nrf2)$geneId)
pathway_name_SWitA_Nrf2<- setReadable(pathway_SWitA_Nrf2, 'org.Hs.eg.db', 'ENTREZID') 


#crear gráficos de las vías enriquecidas de manera individual.
color.params = list(foldChange = as.data.frame(peakAnno_SWitA_p65)$geneId, edge = T)
cex.params = list(category_label = 0.9, gene_label = 0.65)

dotplot(pathway_CN_p65, title = "Condici�n CN Anticuerpo p65")
barplot(pathway_CN_p65, title = "Condici�n CN Anticuerpo p65")
cnetplot(pathway_name_CN_p65, title = "Condici�n CN Anticuerpo p65", circular = T, cex.params= cex.params, color.params = color.params)  + guides(edge_color = "none")

dotplot(pathway_CP_p65, title = "Condici�n CP Anticuerpo p65")
barplot(pathway_CP_p65, title = "Condici�n CP Anticuerpo p65")
cnetplot(pathway_name_CP_p65, title = "Condici�n CP Anticuerpo p65", circular = T, cex.params= cex.params, color.params = color.params)  + guides(edge_color = "none")

dotplot(pathway_Resv_p65, title = "Condici�n Resv Anticuerpo p65")
barplot(pathway_Resv_p65, title = "Condici�n Resv Anticuerpo p65")
cnetplot(pathway_name_Resv_p65, title = "Condici�n Resv Anticuerpo p65",  circular = T, cex.params= cex.params, color.params = color.params)  + guides(edge_color = "none")

dotplot(pathway_SWitA_p65, title = "Condici�n SWitA Anticuerpo p65")
barplot(pathway_SWitA_p65, title = "Condici�n SWitA Anticuerpo p65")
cnetplot(pathway_name_SWitA_p65, title = "Condici�n SWitA Anticuerpo p65",  circular = T, cex.params= cex.params, color.params = color.params)  + guides(edge_color = "none")

dotplot(pathway_CP_Nrf2, title = "Condici�n CP Anticuerpo Nrf2")
dotplot(pathway_Resv_Nrf2, title = "Condici�n Resv Anticuerpo Nrf2")
dotplot(pathway_SWitA_Nrf2, title = "Condici�n SWitA Anticuerpo Nrf2",  circular = T, cex.params= cex.params, color.params = color.params)  + guides(edge_color = "none")


#Functional profiles comparison
peakAnnoList_p65_resv <- lapply(file_p65_resv, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
peakAnnoList_p65_SWitA <- lapply(file_p65_SWitA, annotatePeak, TxDb=txdb,
                                tssRegion=c(-3000, 3000), verbose=FALSE)
peakAnnoList_p65 <- lapply(file_p65, annotatePeak, TxDb=txdb,
                                 tssRegion=c(-3000, 3000), verbose=FALSE)

peakAnnoList_Nrf2 <- lapply(file_Nrf2, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList_p65)
plotAnnoBar(peakAnnoList_Nrf2)

plotDistToTSS(peakAnnoList_p65)
plotDistToTSS(peakAnnoList_Nrf2)


#Recuperar genes anotados para anotacion con Reactome en compareCluster.

genes_p65_resv = lapply(peakAnnoList_p65_resv, function(i) as.data.frame(i)$geneId)
genes_p65_SWitA = lapply(peakAnnoList_p65_SWitA, function(i) as.data.frame(i)$geneId)
genes_Nrf2 = lapply(peakAnnoList_Nrf2, function(i) as.data.frame(i)$geneId)
names(genes_p65_resv) = sub("_", "\n", names(genes_p65_resv))
names(genes_p65_SWitA) = sub("_", "\n", names(genes_p65_SWitA))
names(genes_Nrf2) = sub("_", "\n", names(genes_Nrf2))



compPathway_p65_resv <- compareCluster(geneCluster   = genes_p65_resv,
                           fun           = "enrichPathway",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           readable = T)

compPathway_p65_SWitA <- compareCluster(geneCluster   = genes_p65_SWitA,
                                       fun           = "enrichPathway",
                                       pvalueCutoff  = 0.05,
                                       pAdjustMethod = "BH",
                                       readable = T)


compPathway_Nrf2 <- compareCluster(geneCluster   = genes_Nrf2,
                               fun           = "enrichPathway",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH")

#Graficación de pathways por anticuerpo.
dotplot(compPathway_p65, showCategory = 15, title = "Reactome Pathway Enrichment Analysis")
compPathway_p65_names <- setReadable(compPathway_p65, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(compPathway_p65_names)
viewPathway()

dotplot(compPathway_Nrf2, showCategory = 15, title = "Reactome Pathway Enrichment Analysis")
cnetplot(compPathway_Nrf2)

#Overlap of peaks and annotated genes with Venn diagrams.

genes_venn_p65_resv= lapply(peakAnnoList_p65_resv, function(i) as.data.frame(i)$geneId)
genes_venn_p65_SWitA= lapply(peakAnnoList_p65_SWitA, function(i) as.data.frame(i)$geneId)
genes_venn_p65= lapply(peakAnnoList_p65, function(i) as.data.frame(i)$geneId)

genes_venn_Nrf2= lapply(peakAnnoList_Nrf2, function(i) as.data.frame(i)$geneId)

vennplot(genes_venn_p65, by= "ggVennDiagram") + scale_fill_gradient(low = "white", high = "red")
vennplot(genes_venn_p65_resv, by= "ggVennDiagram") + scale_fill_gradient(low = "white", high = "red")
vennplot(genes_venn_p65_SWitA, by= "ggVennDiagram") + scale_fill_gradient(low = "white", high = "red")
vennplot(genes_venn_Nrf2, by= "ggVennDiagram") + scale_fill_gradient(low = "white", high = "red")

intersect_p65 <- process_region_data(Venn(genes_venn_p65))
intersect_Nrf2 <- process_region_data(Venn(genes_venn_Nrf2))



#Recuperación de las listas de genes por intersección. 

file_conn <- file("gene_name_list_Nrf2.csv", "w") #Cambiar código para cada dataset
for (i in 1:nrow(intersect_Nrf2)) {
  gene_list <- intersect_Nrf2$item[[i]]  # Obtener la lista de genes de la fila actual
  gene_names <- mapIds(org.Hs.eg.db, keys = gene_list, column = "SYMBOL", keytype = "ENTREZID")
  gene_list_str <- paste(gene_names, collapse = ", ")
  writeLines(paste(i, gene_list_str), file_conn)
}
close(file_conn)

#Escribir csv de los pathways. Se debe cambiar para p65.
write.table(as.data.frame(compPathway_Nrf2@compareClusterResult), file="pathway_genes_enrichment_Nrf2.csv", quote = F, row.names = F, sep = ",")

#Recuperar los genes enriquecidos en cada pathway. Se debe cambiar el dataset para p65. 
pathway_analisis <- as.data.frame(compPathway_Nrf2@compareClusterResult)
genes_pathway_analisis <- pathway_analisis[,c(1,9)]
genes_entrez_split <- strsplit(genes_pathway_analisis$geneID, "/")
genenames_list_pathway_Nrf2<- lapply(genes_entrez_split, function(genes) {
  mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENTREZID")
})
genenames_list_pathway_Nrf2 <- sapply(genenames_list_pathway_Nrf2, function(x) paste(x, collapse = "/"))
write.table(as.data.frame(genenames_list_pathway_Nrf2), file = "genenames_pathway_Nrf2.csv", row.names = FALSE, sep = ",")

#Crear tsv de los genes, nombres y símbolos.
write.table(as.data.frame(as.data.frame(peakAnno_CN_p65)[c("SYMBOL", "geneId", "GENENAME")]), file="peakAnno_CN_p65.tsv", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(as.data.frame(peakAnno_CN_Nrf2)[c("SYMBOL", "geneId", "GENENAME")]), file="peakAnno_CN_Nrf2.tsv", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(as.data.frame(peakAnno_CP_p65)[c("SYMBOL", "geneId", "GENENAME")]), file="peakAnno_CP_p65.tsv", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(as.data.frame(peakAnno_CP_Nrf2)[c("SYMBOL", "geneId", "GENENAME")]), file="peakAnno_CP_Nrf2.tsv", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(as.data.frame(peakAnno_Resv_p65)[c("SYMBOL", "geneId", "GENENAME")]), file="peakAnno_Resv_p65.tsv", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(as.data.frame(peakAnno_Resv_Nrf2)[c("SYMBOL", "geneId", "GENENAME")]), file="peakAnno_Resv_Nrf2.tsv", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(as.data.frame(peakAnno_SWitA_p65)[c("SYMBOL", "geneId", "GENENAME")]), file="peakAnno_SWitA_p65.tsv", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(as.data.frame(peakAnno_SWitA_Nrf2)[c("SYMBOL", "geneId", "GENENAME")]), file="peakAnno_SWitA_Nrf2.tsv", quote = F, row.names = F, sep = "\t")

#Crear df de los genes, nombres y símbolos.
CN_p65_genes <- as.data.frame(as.data.frame(peakAnno_CN_p65)[c("geneId", "SYMBOL", "GENENAME")])
CN_Nrf2_genes <- as.data.frame(as.data.frame(peakAnno_CN_Nrf2)[c("geneId", "SYMBOL", "GENENAME")])
CP_p65_genes <- as.data.frame(as.data.frame(peakAnno_CP_p65)[c("geneId", "SYMBOL", "GENENAME")])
CP_Nrf2_genes <- as.data.frame(as.data.frame(peakAnno_CP_Nrf2)[c("geneId", "SYMBOL", "GENENAME")])
Resv_p65_genes <- as.data.frame(as.data.frame(peakAnno_Resv_p65)[c("geneId", "SYMBOL", "GENENAME")])
Resv_Nrf2_genes <- as.data.frame(as.data.frame(peakAnno_Resv_Nrf2)[c("geneId", "SYMBOL", "GENENAME")])
SWitA_p65_genes <- as.data.frame(as.data.frame(peakAnno_SWitA_p65)[c("geneId", "SYMBOL", "GENENAME")])
SWitA_Nrf2_genes <- as.data.frame(as.data.frame(peakAnno_SWitA_Nrf2)[c("geneId", "SYMBOL", "GENENAME")])


###########################################################################
#ChIP peak data set comparison

promoter2 <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList_p65 <- lapply(file_p65, getTagMatrix, windows=promoter2)
plotAvgProf(tagMatrixList_p65, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList_p65, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")

#normal method
plotPeakProf2(file_p65, upstream = 3000, downstream = 3000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb,
              facet = "row")

## binning method 
plotPeakProf2(file_p65, upstream = 3000, downstream = 3000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb,
              facet = "row", nbin = 800)

#Profile of several ChIP peak data binding to body region
plotPeakProf2(file_p65, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)


