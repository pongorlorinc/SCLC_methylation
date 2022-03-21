library(rtracklayer)
library(GenomicRanges)
library(data.table)

# This function scores the average of bw signal in GRanges coordinates defined
# in the "coords" variable. Returns vector of scores with same length as "coords"
ScoreGrangesBWmeanMean = function(coords, bw) {
  overs = as.data.frame(findOverlaps(coords, bw))
  overs$score = bw[overs$subjectHits]$score
  overs.dt = data.table(overs)
  oversum = overs.dt[,list(score = mean(score)), by='queryHits']
  meanvals = rep(NA, length(coords))
  meanvals[oversum$queryHits] = oversum$score
  return (meanvals)
}

# subtract features
subtractOverlap = function(coord, bed) {
  overs = findOverlaps(query = coord, subject = bed)
  return(coord[-unique(queryHits(overs))])
}

bigwig_sample_sheet = "Cell_line_bigwigs.tsv" # sample sheet
annotation_file = "../data/UCSC-HG19-RefSeq_Genes-refGene-2019-07-10.txt" # gene annotation
CpG_island_bed = "../data/hg19.cpg.bed" # CpG island bed file
promoter_extend = 2500

# row names are genes, column names are cell lines.
# Cell line names have to match names for cells for methylation data
# The expression file used was downloaded from the SCLC-CellMinerCDB website
gene_expression_file = "../data/data_nciSclc_exp.txt"

bigwig_samples = read.table(bigwig_sample_sheet, sep = "\t", header = T)
rownames(bigwig_samples) = bigwig_samples[,1]

# import annotation file
annotation = read.table(annotation_file, sep = "\t", header = T)
genes.df = annotation[,c("chrom", "txStart", "txEnd", "strand", "name2", "name")]
colnames(genes.df) = c("chrom", "txStart", "txEnd", "strand", "gene", "tx_name")
genes.df$long_name = paste0(genes.df$gene, "_", genes.df$tx_name, "_", genes.df$chrom, ":", genes.df$txStart, "-", genes.df$txEnd)
genes = makeGRangesFromDataFrame(genes.df, keep.extra.columns = T, seqnames.field = "chrom", start.field = "txStart", end.field = "txEnd", strand.field = "strand")

# Get gene promoters (all transcripts)
gene_promoters = promoters(genes, upstream = promoter_extend, downstream = promoter_extend)
# Import CpG island bed file
cpgs = rtracklayer::import(CpG_island_bed)


# Import bigwig files, store each in "bws" list
bws = list()
for(i in 1:nrow(bigwig_samples)) {
  print(paste0("Importing ", i, " of ", nrow(bigwig_samples)))
  bws[[i]] =   import(bigwig_samples$bigwig[i])
}

quantified.df = data.frame(matrix(ncol = nrow(bigwig_samples), nrow = length(genes)))
colnames(quantified.df) = bigwig_samples$sample
rownames(quantified.df) = genes$long_name
for(i in 1:length(innames)) {
  bw = bws[[i]]
  bw = subtractOverlap(bw, gene_promoters) # remove promoter probes
  bw = subtractOverlap(bw, cpgs) # remove CpG island probes
  quantified.df[,i] = ScoreGrangesBWmeanMean(genes, bw) # score bigwigs
}

write.table(quantified.df, "Quantified_gene_body_methylation_all_transcripts.tsv", col.names = NA, row.names = T, sep = "\t", quote = F)
write.table(genes.df, "Full_gene_annotation.tsv", col.names = NA, row.names = T, sep = "\t", quote = F)

# IF expression file exists, then calculate correlation of gene transcript methylation with gene expression
if(file.exists(gene_expression_file)) {
  # Import expression file
  expdat = read.table(gene_expression_file, header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
  
  # remove rows and columns with NA only
  expdat = expdat[,colSums(is.na(expdat))<nrow(expdat)]
  expdat = expdat[rowSums(is.na(expdat))<ncol(expdat),]
  
  # Get common cell names
  common_cells = intersect(colnames(quantified.df), colnames(expdat))
  
  # Get common genes
  common_genes = intersect(genes$gene, rownames(expdat))
  genes.df_selected = genes.df[genes.df$gene %in% common_genes,]
  genes.df_selected$expression_correlation = NA
  
  # calculate correlation of transcript methylation with gene expression
  genes.df_selected$expression_correlation = sapply(seq(1:nrow(genes.df_selected)), function(i) {
    if(i %% 5000 == 0) {
      print(i)
    }
    cor(as.numeric(quantified.df[genes.df_selected$long_name[i],common_cells]), as.numeric(expdat[genes.df_selected$gene[i], common_cells]), use = "pairwise.complete.obs")
  })
  
  # sort genes based on correlation
  genes.df_selected = genes.df_selected[order(genes.df_selected$expression_correlation, decreasing = T),]
  
  # get best transcript by removing duplicate gene names
  genes.df_selected.best_hits = genes.df_selected[!duplicated(genes.df_selected$gene),]
  
  quantified.df_assigned = quantified.df[genes.df_selected.best_hits$long_name,]
  rownames(quantified.df_assigned) = genes.df_selected.best_hits$gene
  
  write.table(quantified.df_assigned, "Quantified_gene_body_methylation_best_transcript_selected.tsv", col.names = NA, row.names = T, sep = "\t", quote = F)
  write.table(genes.df_selected.best_hits, "Selected_transcript_gene_annotation.tsv", col.names = NA, row.names = T, sep = "\t", quote = F)
}


#gene = "YAP1"
#df = data.frame(body = as.numeric(quantified.df_assigned[gene,common_cells]),
#                exp = as.numeric(expdat[gene,common_cells]))

#library(ggpubr)

#ggplot(df, aes(x = body, y = exp)) +
#  geom_point( size = 3) +
#  xlim(0,1) +
#  stat_cor(method="pearson") +
#  geom_smooth(method=lm, se=FALSE, color = "grey") +
#  ylab("Expression") +
#  xlab("Body methylation") +
#  theme_bw() +
#  theme(aspect.ratio=3/3)
