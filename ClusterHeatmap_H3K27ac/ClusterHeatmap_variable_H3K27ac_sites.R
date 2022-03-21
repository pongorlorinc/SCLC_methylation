library(ggplot2)
library(umap)
library(ComplexHeatmap)
library(rtracklayer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(sva)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

getGeneAnnotation_from_txdb = function(txdatabase = txdb, orgdb = "org.Hs.eg.db") {
  geneann = genes(txdatabase)
  mappings = bitr(geneann$gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb=orgdb)
  geneann = geneann[geneann$gene_id %in% mappings$ENTREZID]
  names(geneann) = geneann$gene_id
  geneann = geneann[mappings$ENTREZID]
  geneann$SYMBOL = mappings$SYMBOL
  promoterann = promoters(geneann, upstream = 2500, downstream = 2500)
  promoterann$SYMBOL  = geneann$SYMBOL
  return(promoterann)
}

norm_counts = read.table("Normalized_variable_enhancer_sites.txt", header = T, row.names = 1, sep = "\t", check.names = F)

# NAPY colors
NAPY_colors = c("ASCL1" = "#FF3100", "NEUROD1" = "#4472C4", "POU2F3" = "#91D050", "YAP1" = "#e0e026")

# get cell annotation
cell_annotation = read.table("Cell_line_annotation.tsv", header = T, sep = "\t", check.names = F)
rownames(cell_annotation) = cell_annotation$global
cell_annotation = cell_annotation[colnames(norm_counts),]

kmthist = kmeans(dist(norm_counts), centers = 4)


ha = HeatmapAnnotation(NAPY = cell_annotation$NAPY, col = list(NAPY = NAPY_colors))
Heatmap(norm_counts, 
        show_row_names = F, 
        top_annotation = ha, 
        column_split = cell_annotation$NAPY,
        row_split = kmthist$cluster,
        use_raster = T) 

