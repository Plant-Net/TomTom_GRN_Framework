#############################
#   DEA for tomtom project  #
#############################

library("DESeq2")
library("ggplot2")
library("stringr")
library("tidyr")
library("dplyr")
library("gplots")
library("RColorBrewer")
library("edgeR")
library("IHW")
source("DEA/Scripts.R")
setwd("DEA")

# Load data
data <- read.delim(file = "../Data/tomato_raw_feature_counts_20240912.tsv", row.names="Geneid", sep='\t')
data <- data %>%
  mutate(across(matches("PMPS23"), .names = "{.col}_copy")) %>%
  rename_with(~ str_replace_all(., "PMPS23", "PM"), matches("_copy$")) %>%
  rename_with(~ str_replace_all(., "PMPS23", "PS23"), matches("PMPS23")) %>%
  rename_with(~ str_replace(., "_copy$", ""), matches("_copy$"))

# Extract metadata
col_names <- colnames(data)
split_names <- str_split_fixed(col_names, "_", 5)
metadata <- data.frame(
  species = split_names[, 1],
  type = split_names[, 2],
  exposition = split_names[, 3],
  location = split_names[, 4],
  number = split_names[, 5]
)

# Define species list
species_list <- c("BC", "CF", "PI", "PM", "PS23", "MI_168hpi", "MI_336hpi")  # Add all species of interest here

# Loop through each species and perform DEA
for (species in species_list) {
  cat("----Processing species:", species, "\n")
  # Check if the species has time points (e.g., based on the 'exposition' column)
  if (str_detect(species,"hpi")) {
    cat("Species has time points, filtering data for time points...\n")
    hpi <- str_split(species, "_")[[1]][2]  # Extract the time point (e.g., "168hpi" or "336hpi")
    species_alone <- str_split(species, "_")[[1]][1]  # Extract the base species name (e.g., "MI")

    #retrieve the data for the species
    res <- create_split_data(data, species_alone)
    data_res <- res$selected_data
    metadat_res <- res$selected_metadata
    data_long_res <- res$selected_data_long
    # Use the function for species with time points
    data_res_hpi <- data_res %>% select(matches(hpi))
    metadat_res_hpi <- metadat_res %>% filter(exposition == hpi)
    data_long_res_hpi <- data_res_hpi %>%
      mutate(Gene = rownames(data_res_hpi)) %>%
      relocate(Gene, .before = 1) %>%
      pivot_longer(cols = -Gene,  # Keep 'Gene' column, pivot all others
                   names_to = "Sample", 
                   values_to = "Expression")
    
    # Replace the original data with the filtered data
    data_res <- data_res_hpi
    metadat_res <- metadat_res_hpi
    data_long_res <- data_long_res_hpi
  } else {
    # Filter data for the current species
    res <- create_split_data(data, species)
    data_res <- res$selected_data
    metadat_res <- res$selected_metadata
    data_long_res <- res$selected_data_long
  }
  
  # Create a directory for the species without "timepoints" in the name
  species_dir <- file.path(getwd(), species)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(species_dir)) {
    dir.create(species_dir)
  }

  # ## DEA
  print("Creating DESeqDataSet...")
  dds <- DESeqDataSetFromMatrix(countData = data_res, colData = metadat_res, design = ~type)
  dds <- estimateSizeFactors(dds)
  
  # ##exploratory plots
  print("Creating rlog transform...")
  rld <- rlog(dds)

  rlog_long <- as.data.frame(assay(rld)) %>%
  mutate(Gene = rownames(rld)) %>%
  relocate(Gene, .before = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression")
  # Plot rlog transformed data
  p <- plot_rnaseq_distributions(rlog_long)
  ggsave(paste0(species_dir,"/distrib_raw_data_rlog_norm.png"), plot = p, width = 15, height = 8, dpi = 500)

  #samples distance
  sampleDists <- dist( t( assay(rld) ) )  
  sampleDistMatrix <- as.matrix( sampleDists ) 
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
  hc <- hclust(sampleDists) 
  pdf(file = paste0(species_dir, "/sample_distance.pdf"), width = 20, height = 10);
  heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=colors, margins=c(2,10),labCol=FALSE )
  dev.off()

  ## PCA plot using ggplot2
  data_pca <- plotPCA(rld, intgroup = "type", returnData = TRUE)
  percentVar <- round(100 * attr(data_pca, "percentVar"))

  # Create PCA plot with ggplot2
  pca_plot <- ggplot(data_pca, aes(x = PC1, y = PC2, color = type, label = name)) +
    geom_point(size = 4) +
    geom_text(hjust = 0.5, vjust = -0.5, size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(legend.position = "right")

  # Save the PCA plot
  ggsave(filename = paste0(species_dir, "/pca_raw_data_vst_norm.png"), plot = pca_plot, width = 12, height = 9, dpi = 300)

  #extract deseq results
  dds <- DESeq(dds)


  res_IHW_05 <- results(dds, filterFun=ihw, alpha=0.05) #facultative, use IHW correction

  write.table(res_IHW_05, paste0(species_dir, "/DEA_Results_ihw.tsv"), sep = "\t")


  # MA plot

  # DESeq2::plotMA(res_IHW_05, ylim=c(-5,5))
  # abline(h = 1,col='red',lwd=4)
  # text(-10,1,'1',col='red')
  # abline(h = -1,col='red',lwd=4)
  # text(-10,-1,'-1',col='red')

  # # #plot dispersion

  # plotDispEsts(dds)


  # # #pvalues plot

  # hist(res_IHW_05$pvalue[res_IHW_05$baseMean > 1], breaks=20, col="grey50", border="white")

  # hist(res_IHW_05$log2FoldChange[res_IHW_05$baseMean > 1], breaks=20, col="grey50", border="white")
  # abline(v=1,col='red')
  # abline(v=-1,col='red')

  # # ##Get the results and do some plots

  # # #heatmap of top DE genes profiles
  # topVarGenes <- head(order(-rowVars(assay(rld))),35)
  # colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255) 
  # sidecols <- c("grey","dodgerblue")[ rld$type ] 
  # mat <- assay(rld)[ topVarGenes, ] 
  # mat <- mat - rowMeans(mat) 

  # heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols, labRow=FALSE, mar=c(10,2), scale="row",cexCol  = 2)

  # # ###
  # topGene <- rownames(res_IHW_05)[which.min(res_IHW_05$padj)]
  # col_dds<-colData(dds)
  # data <- plotCounts(dds, gene=topGene, intgroup=c("type"), returnData=TRUE)

  # ggplot(data, aes(x=type, y=count, fill=type)) + scale_y_log10() + geom_dotplot(binaxis="y", stackdir="center")

}
