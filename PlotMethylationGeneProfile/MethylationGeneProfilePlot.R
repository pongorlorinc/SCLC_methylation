library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)

ScoreGrangesBWmeanMean = function(coords, bw) {
  overs = as.data.frame(findOverlaps(coords, bw))
  overs$score = bw[overs$subjectHits]$score
  nhits =  table(overs$queryHits)
  overs.dt = data.table(overs)
  oversum = overs.dt[,list(score = sum(score)), by='queryHits']
  oversum$score = oversum$score / as.numeric(nhits[as.character(oversum$queryHits)])
  meanvals = rep(NA, length(coords))
  meanvals[oversum$queryHits] = oversum$score
  return (meanvals)
}

# This function is used to bin the coordinates
# intput is:
#   bw: methylation bigwig in GRanges
#   coord: coordinates in GRanges
#   nbins: number of bins
BinCoordinates = function(bw, coord, nbins, remove_NA = T) {
  coord = coord[width(coord) > nbins]
  coord.tiles = unlist(tile(coord, n = nbins))
  coord.tiles$bin = rep(1:nbins, length(coord))
  coord.tiles$name = rep(1:length(coord), each = nbins)
  names(coord.tiles) = NULL
  coord.tiles$score = ScoreGrangesBWmeanMean(coords = coord.tiles, bw = methbw)
  
  if(remove_NA) {
    coord.tiles = coord.tiles[!is.na(coord.tiles$score)]
  }
  
  return (as.data.frame(coord.tiles))
}


# size of upstream and downstream regions
downstream = 2500
upstream = 2500

# RNA expression file, used to separate genes into quantiles
rna_expression_file = "../../data/data_utsw_xsq.txt"

# Gene annotation file
gene_annotation_file = "../../data/UCSC-HG19-RefSeq_Genes-refGene-2019-07-10.txt"

# bigwig file to plot
meth_bigwig_file = "../idat2bw/methylation_bigwigs/NCI-H69.850k_probes.bw"

# CpG island file (used to select transcript)
cpg_file = "../idat2bw/hg19.cpg.bed"

# importing methylation bigwig
methbw = rtracklayer::import(meth_bigwig_file)

# importing cpgs
cpgs = rtracklayer::import(cpg_file)

# importing RNA
rnaseq = read.table(rna_expression_file, header = T, row.names = 1, sep = "\t", stringsAsFactors = F, check.names = F)
rna_means = rowMeans(rnaseq, na.rm = T) # calculate row means
rnadf = data.frame(name = rownames(rnaseq), mean = as.numeric(rna_means)) # create data frame with gene and mean expression
rownames(rnadf) = rownames(rnaseq)
rnadf$quantile = ifelse(rnadf$mean > 0.5, 1, 0) # mark non expressed genes (FPKM <= 1) as 0 quantile
rnadf[rnadf$mean > 1, "quantile"] = ntile(rnadf[rnadf$mean > 1, "mean"], 4) # split expressed genes into quantiles

# import annotation file
annotation = read.table(gene_annotation_file, sep = "\t", header = T)
annotation = annotation[,c("chrom", "txStart", "txEnd", "strand", "name2")]
colnames(annotation) = c("chr", "start", "end", "strand", "gene")
annotation.gr = makeGRangesFromDataFrame(annotation, keep.extra.columns = T)

# getting transcript with closest CpG island
cpgdistance = as.data.frame(distanceToNearest(promoters(annotation.gr, downstream = downstream, upstream = upstream), cpgs))
annotation.gr$cpgdistance = NA
annotation.gr[cpgdistance$queryHits]$cpgdistance = cpgdistance$distance
annotation = as.data.frame(annotation.gr)
annotation = annotation[order(annotation$cpgdistance, decreasing = F),]
annotation = annotation[!is.na(annotation$cpgdistance),]
annotation = annotation[!duplicated(annotation$gene),]
annotation.gr = makeGRangesFromDataFrame(annotation, keep.extra.columns = T)
annotation.gr = sort(annotation.gr)
annotation.gr = annotation.gr[width(annotation.gr) > upstream] # removes short genes

# obtaining gene promoter coordinates
geneproms = promoters(annotation.gr, downstream = downstream, upstream = upstream)

# obtaining gene body coordinates (also removes promoter area from gene body)
genebody = as.data.frame(annotation.gr)
genebody$start = ifelse(genebody$strand == "+", genebody$start + upstream, genebody$start)
genebody$end = ifelse(genebody$strand == "-", genebody$end - upstream, genebody$end)
genebody = makeGRangesFromDataFrame(genebody, keep.extra.columns = T)
genebody$bodylength = width(genebody)
geneproms$length = width(geneproms)

# obtaining gene downstream coordinates
genedownstream = as.data.frame(genebody)
genedownstream$nstart = ifelse(genedownstream$strand == "+", genedownstream$end, genedownstream$start - downstream)
genedownstream$nend = ifelse(genedownstream$strand == "+", genedownstream$end + downstream, genedownstream$start)
genedownstream$start = genedownstream$nstart
genedownstream$end = genedownstream$nend
genedownstream$nstart = NULL
genedownstream$nend = NULL
genedownstream = makeGRangesFromDataFrame(genedownstream, keep.extra.columns = T)

# obtaining gene upstream coordinates
geneupstream = as.data.frame(geneproms)
geneupstream$nstart = ifelse(geneupstream$strand == "+", geneupstream$start - upstream, geneupstream$end)
geneupstream$nend = ifelse(geneupstream$strand == "+", geneupstream$start, geneupstream$end + upstream)
geneupstream$start = geneupstream$nstart
geneupstream$end = geneupstream$nend
geneupstream$nstart = NULL
geneupstream$nend = NULL
geneupstream = makeGRangesFromDataFrame(geneupstream, keep.extra.columns = T)

# gene profile will use these bins
up_and_downstraem_bins = 10 # no. of bins for upstream and downstream area
promoterbins = 25 # no. of bins for promoter area
bodybins = 50 # no of bins for gene body

quants = sort(unique(rnadf$quant)) # get quantiles
plots = list() # store gene plots in a list
for(i in 1:length(quants)) {
  genes = rnadf[rnadf$quant == quants[i], "name"] # get gene name for given quantile
  dfs = list()
  
  # create bin matrix for upstream area
  dfs[[1]] = BinCoordinates(bw = methbw, coord = geneupstream[geneupstream$gene %in% genes], nbins = up_and_downstraem_bins)
  
  # create bin matrix for promoter area
  dfs[[2]] = BinCoordinates(methbw, geneproms[geneproms$gene %in% genes], promoterbins)
  
  # create bin matrix for gene body area
  dfs[[3]] = BinCoordinates(methbw, genebody[genebody$gene %in% genes], bodybins)
  
  # create bin matrix for downstream area
  dfs[[4]] = BinCoordinates(methbw, genedownstream[genedownstream$gene %in% genes], up_and_downstraem_bins)
  
  # offset the bin coordinate to merge the plots
  dfs[[2]]$bin = dfs[[2]]$bin + up_and_downstraem_bins
  dfs[[3]]$bin = dfs[[3]]$bin + up_and_downstraem_bins + promoterbins
  dfs[[4]]$bin = dfs[[4]]$bin + up_and_downstraem_bins + bodybins + promoterbins
  
  dfs[[1]]$type = "upstream"
  dfs[[2]]$type = "promoter"
  dfs[[3]]$type = "gene body"
  dfs[[4]]$type = "downstream"
  
  df = dfs[[1]]
  for(j in 2:4) {
    df = rbind(df, dfs[[j]])
  }
  
  df$quantile = quants[i]
  
  # plot gene profile
  plots[[i]] = ggplot(df, aes(x = factor(bin), y = score, fill = type)) +
    geom_boxplot(outlier.shape = NA, size = .2) +
    ylim(0,1) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank()) +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
}

do.call("grid.arrange", c(plots, ncol=1))

