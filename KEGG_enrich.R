suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("argparse"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("clusterProfiler"))

parse_arg <- function() {
    parser <- ArgumentParser(description = "KEGG enrichment analysis by searching KEGG database under the name of 'zma'(Zea mays). ")
    parser$add_argument("-gene", "--gene", help = "the differential expressed genes from DESeq2 output. ")
    parser$add_argument("-geneID", "--geneID", help = "table with gene label matched to gene ID. ")
    parser$add_argument("-kegg2name", "--kegg2name", help = "KEGG term to name (description) table. ")
    parser$add_argument("-kegg2gene", "--kegg2gene", help = "KEGG term to gene id table (one line one id). ")
    parser$add_argument("-count", "--count", help = "read count. ") # reconstruct the gene universe background
    parser$add_argument("-O", "--output", default = "default", help = "the output name")
    parser$parse_args()
}

###
argv <- parse_arg()
gene <- argv$gene
ID <- argv$geneID # kegg enrichment only support ncbi-id
output <- argv$output
count <- argv$count
kegg2name <- argv$kegg2name
kegg2gene <- argv$kegg2gene
###
if (output == "default") {
    pre <- strsplit(basename(gene), ".txt")[[1]]
    pre <- paste0(pre, ".KEGG")
    output <- pre
} else {
    output <- output
}
print(output)

gene.df <- read.table(gene, sep = "\t", header = T, row.names = 1)
ID.df <- read.table(ID, sep = "\t", header = T, row.names = 1, fill = T) # gene ID and gene label transformation
interest.ID <- ID.df[row.names(gene.df), "Gene_ID", drop = F]
interest.id <- interest.ID[, "Gene_ID", drop = T]
clean.interest.id <- interest.id[!is.na(interest.id)]
### universal background
count <- read.table(count, sep = "\t", header = T, row.names = 1)
out.gene <- rownames(count[rowMeans(count < 10) == 1, ])
out.id <- ID.df[out.gene, "Gene_ID"]
rem.uni <- length(out.id)
print(paste0("Remove ", rem.uni, " gene ID as low expression. "))
#
kegg2gid <- read.table(kegg2gene, sep = "\t", header = T)
out.row <- which(kegg2gid$Gene_ID %in% out.id)
kegg2gid <- kegg2gid[-out.row, ]
#
kegg2nm <- read.table(kegg2name, sep = "\t", header = T, quote = "")

###
 # ID for DEGs, some DEGs don't have ID (NA)
# write.table(out.ID, file = paste0(out, "_ID.txt"), sep = "\t", row.names = T, quote = F)
naid <- length(interest.id[is.na(interest.id)])
print(paste0("NA Value for DEGs ", naid))
###
enrich <- enricher(clean.interest.id, TERM2GENE = kegg2gid, TERM2NAME = kegg2nm, pvalueCutoff = 0.01, pAdjustMethod = "BH")
ekegg.df <- as.data.frame(enrich)
ekegg.df <- ekegg.df[ekegg.df$p.adjust < 0.01, ]

if (nrow(ekegg.df) > 0) {
    print(paste0(nrow(ekegg.df), " enriched for this gene set. "))
    write.table(ekegg.df, file = paste0(output, ".txt"), sep = "\t", row.names = F, quote = F)
} else {
    print("No enrichment for this gene set. ")
}

# NOTE: the downloaded kegg id is not complete and not the lasted kegg pathways.
# when performing enrichKEGG, and kegg pathways are the most updated.