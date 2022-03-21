library(rtracklayer)
library(readxl)
library(edgeR)
library(sva)

promoter_upstreams = c(2500, 5000, 10000, 15000)
promoter_downstreams = c(2500, 5000, 10000, 15000)

# Raw coverage file
raw_count_matrix_histone_file = "data/raw_coverages.tsv"

# table that maps sample names from the raw_coverages.tsv to cell line name, and Z-scored SCLC-CellminerCDB name
# Last column has the batch variable
cell_line_name_file = "data/Name_pairs.xlsx"

# gene annotation file
gene_annotation_file = "../data/UCSC-HG19-RefSeq_Genes-refGene-2019-07-10.txt"

# Genes selected with highes SD (standard deviation) in expression based on RNAseq data from SCLC-CellMinerCDB
high_sd_genes = read.table("data/High_SD_genes_UTSW_xsq.tsv", sep = "\t")[,1]

# Global Z-scored expression file for genes from SCLC-CellMinerCDB
expdata = read.table("../data/data_genSclc_exp.txt", sep = "\t", row.names = 1, header = T, check.names = F)


# Import raw counts
raw_counts = read.table(raw_count_matrix_histone_file, 
                        row.names = 1, 
                        header = T, 
                        check.names = F)

# Import names of samples
#name	\t cell_line \t global
#H1048_SCLC_H3K27ac_ChIP_0917.sorted.dedup.bam \t NCI-H1048 \t NCI-H1048
#H1092_SCLC_H3K27ac_ChIP_0917.sorted.dedup.bam \t NCI-H1092 \t NCI-H1092
cellpairs = as.data.frame(read_xlsx(cell_line_name_file))

# Import gene annotation
annotation_full = read.table(gene_annotation_file, header = T, sep = "\t")
annot.df = annotation_full[,c("chrom", "txStart", "txEnd", "strand", "name2", "name")]
colnames(annot.df) = c("chr", "start", "end", "strand","gene", "transcript")
annot = makeGRangesFromDataFrame(annot.df, keep.extra.columns = T)
names(annot) = seq(1,length(annot))

# Order columns to match name pairing
rownames(cellpairs) = cellpairs$name
cellpairs = cellpairs[colnames(raw_counts),]
cellpairs = cellpairs[!is.na(cellpairs$name),]
raw_counts = raw_counts[,cellpairs$name]

# Get coordinate lengths for TMM normalization
peak_coords = GRanges(rownames(raw_counts))
gene_lengths = data.frame(GeneID = rownames(raw_counts), Length = width(peak_coords))
rownames(gene_lengths) = gene_lengths$GeneID

# Perform normalization
adjusted_counts <- ComBat_seq(as.matrix(raw_counts), batch=cellpairs$batch)
gene_lengths = gene_lengths[rownames(adjusted_counts),]
GeneDF_EdgeR <- edgeR::DGEList(counts = adjusted_counts, genes = gene_lengths)
GeneDF_Norm  <- edgeR::calcNormFactors(GeneDF_EdgeR, method = 'TMM', )
histone_norm_mat <- as.data.frame(edgeR::rpkm(GeneDF_Norm, normalized.lib.sizes = TRUE, log = FALSE))
histone_norm_mat_log2 = log2(histone_norm_mat + 1)
peak_histone_mean_signal = rowMeans(histone_norm_mat_log2, na.rm = F)


overlap_dfs = list()
geneann = annot
for(i in 1:length(promoter_upstreams)) {
  proms = promoters(x = geneann, upstream = promoter_upstreams[i], downstream = promoter_downstreams[i])
  overs = findOverlaps(query = peak_coords, subject = proms)
  overlapping_genes = proms[unique(subjectHits(overs))]$gene
  overlap_dfs[[i]] = as.data.frame(overs)
  colnames(overlap_dfs[[i]]) = c("peak", "feature")
  overlap_dfs[[i]]$feature = names(proms[overlap_dfs[[i]]$feature])
  overlap_dfs[[i]]$level = i
  overlap_dfs[[i]]$mean = peak_histone_mean_signal[overlap_dfs[[i]]$peak]
}

overlap_df = data.frame()

for(i in 1:length(promoter_upstreams)) {
  if(nrow(overlap_df) == 0) {
    overlap_df = overlap_dfs[[i]]
  } else {
    overlap_df = rbind(overlap_df, overlap_dfs[[i]])
  }
}

overlap_df = overlap_df[order(overlap_df$level, -overlap_df$mean, decreasing = F),]
overlap_df = overlap_df[!duplicated(overlap_df$feature),]
overlap_df$peak_name = names(peak_histone_mean_signal[overlap_df$peak])
overlap_df$gene_name = proms[overlap_df$feature]$gene
overlap_df = overlap_df[!duplicated(overlap_df$gene_name),]
write.table(overlap_df, "Promoter_gene_coordinate_pairs.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

peak_histone_gene_signal = histone_norm_mat_log2[overlap_df$peak_name,]
rownames(peak_histone_gene_signal) = overlap_df$gene_name


cellpairs$cor = NA
for(i in 1:nrow(cellpairs)) {
  common_genes = intersect(rownames(peak_histone_gene_signal), rownames(expdata))
  common_genes = intersect(common_genes, high_sd_genes)
  cellpairs$cor[i] = cor(peak_histone_gene_signal[common_genes, i], expdata[common_genes,cellpairs$global[i]], use = "complete.obs")
}

write.table(cellpairs, "Cell_line_correlation_with_gene_expression_TMM_FPKM.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

selected_cell_pairs = cellpairs
selected_cell_pairs = selected_cell_pairs[order(selected_cell_pairs$cor, decreasing = T),]
selected_cell_pairs = selected_cell_pairs[!duplicated(selected_cell_pairs$cell_line),]
write.table(selected_cell_pairs, "Cell_line_correlation_with_gene_expression_selected_TMM_FPKM.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

peak_histone_gene_signal = peak_histone_gene_signal[,selected_cell_pairs$name]
colnames(peak_histone_gene_signal) = selected_cell_pairs$global
write.table(peak_histone_gene_signal, "H3K27ac_gene_promoter_values_TMM_FPKM.tsv", sep = "\t", quote = F, col.names = NA, row.names = T)
