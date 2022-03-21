library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(umap)

# data frame (table) with cell line names, path to bigwig and NAPY status
# Need 3 columns: "sample", "bigwig" and "NAPY"
#
# sample  bigwig                                                NAPY
# DMS-53	../idat2bw/methylation_bigwigs/DMS-53.850k_probes.bw	NEUROD1
cell_line_bws = read.table("Cell_line_bigwigs_withNAPY.tsv", header = T, sep = "\t")

promoter_extend_bases = 2500 # extension of promoter in bp

# CpG island BED file
cpgs = rtracklayer::import("../data/hg19.cpg.bed") 

# setting up variables for UMAP
custom.config = umap.defaults
custom.config$random_state = 1
custom.config$n_components = 2
custom.config$n_neighbors = 10
custom.config$min_dist = 0.1

# setting row names
rownames(cell_line_bws) = cell_line_bws$sample
# Color of cell lines
napy_colors = c("ASCL1" = "red", "NEUROD1" = "green", "POU2F3" = "blue", "YAP1" = "yellow")

# Importing bigwig files
bws = list()
for(i in 1:nrow(cell_line_bws)) {
  print(paste0("Importing ", i, " of ", nrow(cell_line_bws)))
  bws[[i]] = rtracklayer::import(cell_line_bws$bigwig[i])
}

# Combining bigwig files into a data frame
probe.df = data.frame(matrix(ncol =  nrow(cell_line_bws), nrow = length(bws[[1]])))
refprobe = bws[[1]] # use first bigwig as reference probe set
names(refprobe) = seq(1:length(refprobe)) # name probes
colnames(probe.df) =  cell_line_bws$sample # add column name (defined in sample column)
rownames(probe.df) = names(refprobe) # add probe names

for(i in 1:nrow(cell_line_bws)) {
  overs = findOverlaps(query = bws[[i]], subject = refprobe) # find overlaps
  probe.df[subjectHits(overs),i] = bws[[i]][queryHits(overs)]$score # add overlaps to data frame
}

probe.df = probe.df[complete.cases(probe.df),] # remove rows with NA values
refprobe = refprobe[rownames(probe.df),] # subset reference probe set

# Annotation of probes
refprobe$annotation = "Intergenic"

# getting gene annotating probes
overs = findOverlaps(query = refprobe, subject = genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
refprobe[queryHits(overs)]$annotation = "Gene body"

# getting promoter annotating probes
overs = findOverlaps(query = refprobe, subject = promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = promoter_extend_bases, downstream = promoter_extend_bases))
refprobe[queryHits(overs)]$annotation = "Promoter"

# getting CpG annotating probes
overs = findOverlaps(query = refprobe, subject = cpgs)
refprobe[queryHits(overs)]$annotation = "CpG"

# UMAP of all probes
probe.df.umap = umap(t(probe.df), method = "naive", config = custom.config)
probe.df.umap.plot <- data.frame(x = probe.df.umap$layout[,1], 
                                 y = probe.df.umap$layout[,2], 
                                 name = rownames(probe.df.umap$layout), 
                                 NAPY = cell_line_bws[rownames(probe.df.umap$layout), "NAPY"])

probe.df.umap.plot$color = napy_colors[probe.df.umap.plot$NAPY]
probe.df.umap.plot$color = ifelse(is.na(probe.df.umap.plot$color), "darkgrey", probe.df.umap.plot$color)

pdf("All_probes.umap_with_names.pdf")
ggplot(probe.df.umap.plot, aes(x = x, y = y, label=name)) + 
  geom_point(fill = probe.df.umap.plot$color, size = 4, colour="black", pch=21) +
  geom_text_repel(aes(label=name), size = 2) +
  theme_bw()
dev.off()


# UMAP of CpG probes
probe.df.umap = umap(t(probe.df[refprobe$annotation == "CpG",]), method = "naive", config = custom.config)
probe.df.umap.plot <- data.frame(x = probe.df.umap$layout[,1], 
                                 y = probe.df.umap$layout[,2], 
                                 name = rownames(probe.df.umap$layout), 
                                 NAPY = cell_line_bws[rownames(probe.df.umap$layout), "NAPY"])

probe.df.umap.plot$color = napy_colors[probe.df.umap.plot$NAPY]
probe.df.umap.plot$color = ifelse(is.na(probe.df.umap.plot$color), "darkgrey", probe.df.umap.plot$color)

pdf("CpG_probes.umap_with_names.pdf")
ggplot(probe.df.umap.plot, aes(x = x, y = y, label=name)) + 
  geom_point(fill = probe.df.umap.plot$color, size = 4, colour="black", pch=21) +
  geom_text_repel(aes(label=name), size = 2) +
  theme_bw()
dev.off()

# UMAP of promoter probes
probe.df.umap = umap(t(probe.df[refprobe$annotation == "Promoter",]), method = "naive", config = custom.config)
probe.df.umap.plot <- data.frame(x = probe.df.umap$layout[,1], 
                                 y = probe.df.umap$layout[,2], 
                                 name = rownames(probe.df.umap$layout), 
                                 NAPY = cell_line_bws[rownames(probe.df.umap$layout), "NAPY"])

probe.df.umap.plot$color = napy_colors[probe.df.umap.plot$NAPY]
probe.df.umap.plot$color = ifelse(is.na(probe.df.umap.plot$color), "darkgrey", probe.df.umap.plot$color)

pdf("Promoter_probes.umap_with_names.pdf")
ggplot(probe.df.umap.plot, aes(x = x, y = y, label=name)) + 
  geom_point(fill = probe.df.umap.plot$color, size = 4, colour="black", pch=21) +
  geom_text_repel(aes(label=name), size = 2) +
  theme_bw()
dev.off()

# UMAP of gene body probes
probe.df.umap = umap(t(probe.df[refprobe$annotation == "Gene body",]), method = "naive", config = custom.config)
probe.df.umap.plot <- data.frame(x = probe.df.umap$layout[,1], 
                                 y = probe.df.umap$layout[,2], 
                                 name = rownames(probe.df.umap$layout), 
                                 NAPY = cell_line_bws[rownames(probe.df.umap$layout), "NAPY"])

probe.df.umap.plot$color = napy_colors[probe.df.umap.plot$NAPY]
probe.df.umap.plot$color = ifelse(is.na(probe.df.umap.plot$color), "darkgrey", probe.df.umap.plot$color)

pdf("GeneBody_probes.umap_with_names.pdf")
ggplot(probe.df.umap.plot, aes(x = x, y = y, label=name)) + 
  geom_point(fill = probe.df.umap.plot$color, size = 4, colour="black", pch=21) +
  geom_text_repel(aes(label=name), size = 2) +
  theme_bw()
dev.off()
