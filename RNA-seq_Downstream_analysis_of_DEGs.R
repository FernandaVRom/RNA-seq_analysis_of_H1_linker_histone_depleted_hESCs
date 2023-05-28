##################################################################################
# TITLE - RNA-seq. Downstream analysis of DEGs                                   #
# AUTHOR - Fernanda Vargas-Romero                                                #
# DATE - September 5th, 2022                                                     #
# DESCRIPTION - This R script describes the downstream steps used after obtaining# 
# differentially expressed genes. These steps include the generation of heat maps# 
# and performing functional analyses, such as gene ontology analysis.            #
##################################################################################

############################## IMPORTANT NOTES BEFORE STARTING ###############################
# 1. A "Data" folder has been included into the repository, this will include the processed  #
#    and normalized gene expression data, such as count matrices and metadata, once the work #
#    is ready for publication.                                                               #
# 2. A "Results" folder has been included into the repository. We have used this path to have# 
#    all "write.csv and "save" commands for our different pipelines                          #
##############################################################################################


# HEATMAP GENERATION ----------------------------------------------------------------------------------------
# The next steps are useful to generate heatmaps to visualize the changes in gene expression from different conditions
# (e.g. visualizing changes from SKOs or DKOs)

## Import data -----------------------------------------------------------------------------------------------------------------------------------------------
# You will have to generate your own list of genes to use this command. 
# This list will have to include the ENSEMBL or gene list + the LFC for each of the conditions you are interested in visualizing.
# An example will be included in the data folder once the work is ready for publication.

list1 <-read.csv('./Data/list_of_genes_LFC.csv',
                 header=TRUE, 
                 row.names=1, 
                 check.names=FALSE,
                 )

head(list1)

## Heatmap generation ----------------------------------------------------------------------------------------------------------------------------------------
#load the next libraries

library(pheatmap)
library(stats)

### Generate heatmap----------------------------------------------------------------------------------------------------------------

heatm <- pheatmap(list1, 
                  scale = "row", 
                  cluster_rows = TRUE,
                  color=colorRampPalette(c("blue", "white", "red"))(20),
                  cluster_cols = TRUE,  
                  cutree_cols=3, 
                  cutree_rows=5,
                  show_rownames = FALSE
                  )

### Getting the order of the genes from top to bottom------------------------------------------------------------------------------

order_genes <- rownames(list1)[heatm$tree_row[["order"]]]

# Create a data frame with the order and gene names
order_data <- data.frame(order = 1:length(order_genes), ID = order_genes)

# Save the data frame to a CSV file
write.csv(order_data, './Results/list_of_genes_heatmap_in_order.csv', row.names = FALSE)

### Getting the genes that belong to each cluster --------------------------------------------------------------------------------

# Obtain the cluster assignments

clusters <- cutree(heatm$tree_row, k = 5)

# The k = 5 depends on how many clusters you want to observe on the heatmap. 

# Create a data frame with gene names and cluster assignments
cluster_data <- data.frame(ID = rownames(list1), Cluster = clusters)

# Save the data frame to a CSV file
write.csv(cluster_data, "./Results/list_of_genes_heatmap_clusters.csv", row.names = FALSE)


### Save the heatmap into pdf ----------------------------------------------------------------------------------------------------

save_pheatmap_pdf <- function(x, filename, width=10, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


save_pheatmap_pdf(heatm, './Results/heatmap_clusters.pdf')


# DOWNSTREAM ANALYSIS OF HEATMAP OUTPUT -----------------------------------------------------------------------------------------------------------
# Once you get the file for the order of genes and the genes in each cluster, you can merge both files and get the location of 
# each gene in each cluster along the heatmap

### Obtaining the order of each cluster along the map ----------------------------------------------------------------------------
# Load your files in R in the next order

a1 <- read.csv('./Results/list_of_genes_heatmap_clusters.csv')
head(a1)

a2 <- read.csv('./Results/list_of_genes_heatmap_in_order.csv')
head(a2) 

#merge a1 into a2 and assign them by ID
a3 <- merge(a1, a2, by ="ID")
head(a3)

#save the file 
write.csv(a3, './Results/list_of_genes_heatmap_in_order_and_by_cluster.csv')


### Functional analysis ---------------------------------------------------------------------------------------------------------
# To make a functional analysis, there are different tools. We used both gProfiler or ClusterProfiler. However, the final data
# was obtained from gProfiler (functional analysis for SKO and DKO gene clusters and functional analysis on synergic class on each combination) )

# You will have to generate your own list of genes to use this command. 
# This list will have to include the ENSEMBL + the LFC for each of the conditions you are interested in visualizing.
# An example will be included in the data folder once the work is ready for publication.


#### With gProfiler -------------------------------------------------------------------------------------------------------------
# STEP 1: go to https://biit.cs.ut.ee/gprofiler/
# STEP 2: on g:GOSt upload the list of genes of interest 
# STEP 3: the parameters used for our data are the following
#         In "advanced options" select:
#                                      - All results
#                                      - Only annotated genes
#                                      - Significant threshold Benjamini Hochberg FDR
#                                      - User threshold 0.05
#         In "data sources" select:
#                                      - GO biological process
#STEP 4: after the analysis is done filter the number of terms to the ones between 20-200 or 20-1000


#### With ClusterProfiler -------------------------------------------------------------------------------------------------------

# Load required packages
  
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

# Prepare INPUT 
# Conversion of ENSGID to ENTREZ ID 

# Load your file with LFC and ENSEMBL IDs
data <- read.csv('./Data/list_of_genes_LFC_for_GO_ENSEMBL.csv')

# Convert ENSEMBL IDs to ENTREZ IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = data$ENSGID,
                     column = "ENTREZID",
                     keytype = "ENSEMBL")

# Create a new data frame with LFC and ENTREZ IDs
new_data <- data.frame(ENTREZID = entrez_ids, LFC_Condition_A = data$Condition_A)

# Remove the ENSGID column
new_data$ENSGID <- NULL

# Remove rows with missing or "NA" ENTREZ IDs
new_data <- new_data[!is.na(new_data$ENTREZID), ]

# Reorder the columns (optional)
new_data <- new_data[, c("ENTREZID", "LFC_Condition_A")]

# Write the new file with LFC and ENTREZ IDs
write.csv(new_data, './Results/list_of_genes_ENTREZID_LFC.csv', row.names = FALSE)


# Loading data with ENTREZID + LFC 

Gene_list <- read.csv('./Results/list_of_genes_ENTREZID_LFC.csv', header = TRUE)

# Filter upregulated genes (LFC > 0)
upregulated_genes <- Gene_list[Gene_list$LFC_Condition_A > 0, ]

# Filter downregulated genes (LFC < 0)
downregulated_genes <- Gene_list[Gene_list$LFC_Condition_A < 0, ]

# Perform GO over-representation analysis for upregulated genes
upregulated_GO_results <- enrichGO(gene          = upregulated_genes$ENTREZID,
                                   OrgDb         = org.Hs.eg.db,
                                   keyType       = "ENTREZID",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   minGSSize     = 20, 
                                   maxGSSize     = 200, 
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE
                                   )
# Save GO results for upregulated genes in a .csv file
write.csv(upregulated_GO_results,'./Results/GO_result_upregulated.csv')


# Perform GO over-representation analysis for downregulated genes
downregulated_GO_results <- enrichGO(gene          = downregulated_genes$ENTREZID,
                                   OrgDb         = org.Hs.eg.db,
                                   keyType       = "ENTREZID",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   minGSSize     = 20, 
                                   maxGSSize     = 200, 
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE
                                   )

# Save GO results for downregulated genes in a .csv file
write.csv(downregulated_GO_results,'./Results/GO_result_downregulated.csv')

# Visualize GO enrichment for upregulated genes
Dotp_up <- dotplot(upregulated_GO_results, 
        color = "p.adjust",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 10,
        orderBy = "x",
        label_format = 30,
        )

# Visualize GO enrichment for downregulated genes
Dotp_down <- dotplot(downregulated_GO_results, 
        color = "p.adjust",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 10,
        orderBy = "x",
        label_format = 30,
        )

.# Save the files in pdf
ggsave(filename = "./Results/GO_result_upregulated.pdf", plot = Dotp_up, width = 10, height = 7)
ggsave(filename = "./Results/GO_result_downregulated.pdf", plot = Dotp_down, width = 10, height = 7)

