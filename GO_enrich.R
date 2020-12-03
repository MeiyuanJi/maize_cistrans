suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("dplyr"))

parse_arg <- function() {
    parser <- ArgumentParser(description = "GO enrichment analysis based on hypergeometic distribution. ")
    parser$add_argument("-gene", "--geneList", help = "Interested gene List for this enrichment analysis. ")
    parser$add_argument("-term2name", "--term2name", help = "GO term to name (description) table. ")
    parser$add_argument("-count", "--count", help = "raw count. ")
    parser$add_argument("-term2gene", "--term2gene", help = "GO term to gene id table (one line one id). ")
    parser$add_argument("-O", "--output", default = "default", help = "prefix of output file. ")
    parser$parse_args()
}

###
argv <- parse_arg()
genefile <- argv$geneList
term2name <- argv$term2name
term2gene <- argv$term2gene
output <- argv$output
count <- argv$count

if (output == "default") {
    pre <- strsplit(basename(genefile), ".txt")[[1]]
    pre <- paste0(pre, ".GO")
    output <- pre
} else {
    output <- output
}
print(output)
### remove lowly expressed genes and rebuild universal gene background 
count <- read.table(count, sep = "\t", header = T, row.names = 1)
out.gene <- rownames(count[rowMeans(count < 10) == 1, ])
term2gid <- read.table(term2gene, sep = "\t", header = T)
out.row <- which(term2gid$gene %in% out.gene)
term2gid <- term2gid[-out.row, ]
### enricher() analysis
gene <- read.table(genefile, sep = "\t", header = T, row.names = 1)
#gene.row <- row.names(gene[row.names(gene) %in% out.gene,])
#print(gene.row)
gene <- rownames(gene)
#write.table(gene, file = paste0(pre.gene, ".txt"), sep = "\t", quote = F)
term2nm <- read.table(term2name, sep = "\t", header = T, quote = "")
enrich <- enricher(gene, TERM2GENE = term2gid, TERM2NAME = term2nm, pvalueCutoff = 0.01, pAdjustMethod = "BH")
### 
enrich <- as.data.frame(enrich)
enrich <- enrich[enrich$p.adjust < 0.01, ]
#enrich <- subset(enrich, select = c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count", "geneID"))
if (nrow(enrich) > 0) {
    print(paste0(nrow(enrich), " enriched for this gene set. "))
    write.table(enrich, file = paste0(output, ".txt"), sep = "\t", quote = F, row.names = F)
} else {
    print("No enrichment for this gene set. ")
}


