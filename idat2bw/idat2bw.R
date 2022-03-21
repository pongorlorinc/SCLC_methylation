library(minfi)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)

hg19 = BSgenome.Hsapiens.UCSC.hg19
outdir = "methylation_bigwigs/" # results will go here

input_sample_sheet = "sample_sheet_input.tsv" # sample sheet
# Sample sheet has three coulumns:
# sample	grn_file	rd_file
# NAME    <path>    <path>

if(!dir.exists(outdir)) {
  dir.create(outdir)
}

# import sample sheet
sample_sheet = read.table(input_sample_sheet, header = T, sep = "\t")
# add basename column
sample_sheet$Basename = gsub(pattern = "_Red.idat.gz", replacement = "", x = sample_sheet$rd_file)
# Read data
RGset = read.metharray.exp(targets = sample_sheet)
# Normalize data
normdat = preprocessIllumina(RGset)

# map probes to genome
genome_mapped = mapToGenome(RGset)
rownames(sample_sheet) = basename(sample_sheet$Basename)
# get methylation score
methylation_score = getBeta(genome_mapped)
colnames(methylation_score) = sample_sheet[colnames(methylation_score),1]

# export bigwig files
for(i in 1:ncol(methylation_score)) {
  coords = granges(genome_mapped)
  coords$score = methylation_score[names(coords),i]
  coords = coords[!is.na(coords$score)]
  seqlengths(coords) = seqlengths(hg19)[names(seqlengths(coords))]
  
  outbwfile = paste0(outdir, "/", colnames(methylation_score)[i], ".850k_probes.bw")
  rtracklayer::export(coords, outbwfile)
}