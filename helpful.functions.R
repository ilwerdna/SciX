# These are some functions to help you filter and sort your results

###Find upregulated genes ###
filter.upregulated.genes <- function(data, logfc=0.5, fdr=0.05){
  data <- subset(data, logFC > logfc & FDR < fdr)
  return(data[order(data$logFC, decreasing = TRUE),])
}

### Find downregulated genes ###
filter.downregulated.genes <- function(data, logfc=-0.5, fdr=0.05){
  data <- subset(data, logFC < logfc & FDR < fdr)
  return(data[order(data$logFC, decreasing = FALSE),])
}

### Search for a specific gene ###
find.single.gene <- function(data, target){
  return(subset(data, gene %in% target))
}

### Search for a list of genes ###
find.multiple.genes <- function(data, list){
  return(subset(data, gene %in% list))
}

## Here is an example for working on the output from edgeR which is saved in your analysis folder ###
# Read in file. *Change file path*
lrt <- read.delim('edgeR-LRT.PD.txt')
# Find gene in result file
find.single.gene(lrt, target='IGHM')
# Filter for upregulated genes
upregulated.genes <- filter.upregulated.genes(data=lrt)
# Filter for downregulated genes
downregulated.genes <- filter.downregulated.genes(data=lrt)
# Find a single gene in the upregulated data
find.single.gene(data=upregulated.genes, target='IL1RL1')
# Create a list of genes and search for them in the upregulated data
gene.list <- c('ALB', 'SST', 'CTSV', 'CSF3')
find.multiple.genes(upregulated.genes, list=gene.list)


### Create a volcano plot and label with specific genes ###
volcano.plot.gene.list <- function(data, list, logfc = 0.5, fdr = 0.05) {
  
  threshold <- data$FDR < fdr & abs(data$logFC) > logfc
  data$threshold <- threshold
  
  volcano.plot <- ggplot(data) +
    geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
    geom_text_repel(data=data[data$gene %in% list,], aes(x=logFC, y=-log10(FDR), label=gene), force_pull = 0.1, min.segment.length = 0.05, force=2) +
    ggtitle("Volcano Plot") +
    xlab("log2 fold change") +
    ylab("-log10 FDR") +
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    scale_color_discrete(name = "DEG")
  pdf(paste0('Volcano.Plot.gene.list.pdf'))
  volcano.plot
  dev.off()
  volcano.plot
}
# Example command
volcano.plot.gene.list(lrt, gene.list)                                    

### Boxplot of specific genes ###
# TMM is the normalised data from edgeR that we use for plotting
boxplot.gene.list <- function(tmm, list){
  rownames(tmm) <- gsub('ENSG[0-9]+.[0-9]+\\|', '', rownames(tmm))
  colnames(tmm) <- gsub('[0-9]+', '', colnames(tmm))
  plot.data <- subset(tmm, rownames(tmm) %in% gene.list)
  plot.data <- plot.data[order(rowSums(plot.data), decreasing = T),]
  plot.data <- melt(plot.data)
  colnames(plot.data) <- c('gene', 'condition', 'cpm')
  boxplot.plot <- ggplot(plot.data, aes(x=gene, y=cpm, fill=condition)) + 
    geom_boxplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  pdf(paste0('boxplot.plot.gene.list.pdf'))
  boxplot.plot
  dev.off()
  boxplot.plot
}
# Example command
boxplot.gene.list(tmm, list=gene.list)

### Function to print the pathways identified in GSEA. Result ordered by number of genes in pathway ###
pathways <- function(ego){
  ego.df <- data.frame(ego)
  ego.df[order(ego.df$Count, decreasing = T),1]
}
# Example command
pathways(ego)

### Function to identify genes participating specific GSEA pathway ###
genes.in.pathway <- function(ora, pathway, list=NULL){
  ora.df <- data.frame(ora)
  ora.df <- subset(ora.df, ID %in% pathway)
  genes <- unlist(strsplit(ora.df$geneID, '/'))
  if (!is.null(list)){
    intersected.genes <- intersect(genes, list)
    if (length(intersected.genes) > 0){
      print(genes)
      print('### Specific genes ###')
      return(intersected.genes)
    }
  } else {
    return(genes)
  }
}
# Example command to print genes in an interesting pathway
genes.in.pathway(ego, pathway='REACTOME_NEUTROPHIL_DEGRANULATION')
# Example command to check if specific genes are in pathway
genes.in.pathway(ego, pathway='REACTOME_SIGNALING_BY_INTERLEUKINS', list=c('IL1B', 'CXCL1'))

enrichment.test <- function(data, list, logfc=0.5, fdr=0.05){
  a <- length(which(subset(data, abs(logFC) > logfc & FDR < fdr)$gene %in% list))
  b <- length(which(subset(data, FDR > fdr)$gene %in% list))
  c <- length(which(!(subset(data, abs(logFC) > logfc & FDR < fdr)$gene %in% list)))
  d <- length(which(!(subset(data, FDR > fdr)$gene %in% list)))
  ctable <- matrix(c(a,b,c,d), nrow=2)
  chisq <- chisq.test(ctable)
  return(list(result = chisq, expected.values = chisq$expected, observed.values = chisq$observed, pearson.residuals = chisq$residuals, p.value = chisq$p.value))
}

# Example command to test if X chromosome genes are enriched in upregulated data 
load('../chrX.Rdata')
enrichment.test(lrt, list=rownames(chrX))
