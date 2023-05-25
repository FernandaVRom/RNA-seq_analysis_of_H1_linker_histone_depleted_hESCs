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
                  cutree_rows=2,
                  show_rownames = FALSE
                  )

### Getting the order of the genes from top to bottom------------------------------------------------------------------------------

order_genes <- rownames(list1[heatm$tree_row[["order"]],])

write.csv(order_genes, './Results/list_of_genes_heatmap_in_order.csv')

### Getting the genes that belong to each cluster --------------------------------------------------------------------------------

# The k = 5 depends on how many clusters you want to observe on the heatmap. These should match the optimal number of clusters 
# obtained earlier (see pipeline RNA-seq_quality_assessment_and_processing.Rmd for more)

write.csv(sort(cutree(heatm$tree_row, k=5)), './Results/list_of_genes_heatmap_clusters.csv')

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
# STEP 1 in the .csv files name the columns of genes as "ID"

# STEP 2:load your files in R in the next order

a1 <- read.csv('./Results/list_of_genes_heatmap_clusters.csv')

head(a1)

a2 <- read.csv('./Results/list_of_genes_heatmap_in_order.csv')

head(a2)

#merge a1 into a2 and assign them by ID
a3 <- merge(a1, a2, by ="ID")

head(a3)

#save the file (function needs to be completed)
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


# Load the annotation you will use
# SET THE DESIRED ORGANISM HERE (in this case human)
organism = "org.Hs.eg.db"

library(organism, character.only = TRUE)

# Prepare INPUT 
# Conversion of ENSGID to ENTREZ ID 

keytypes(org.Hs.eg.db) #HELPFUL to see the diffent keys

a1 <- read.delim('./Data/list_of_genes_LFC_for_GO_ENSEMBL.csv', sep = ",")

head(a1)

# You need to specifically tell to keys which column to look into
ENTREZID.list <- mapIds(org.Hs.eg.db, keys = a1$ENSGID, keytype="ENSEMBL", column = "ENTREZID")

head(ENTREZID.list)


# Make a new dataframe so that your ENTREZID goes next to your ENSGID
a2 <- data.frame(a1$ENSGID, ENTREZID.list)

head(a2)

# Save it and transfer the ENTREZ to your DEG table on excel 
write.csv(a2,'./Results/list_of_genes_for_GO_ENSEMBL_to_ENTREZ.csv', row.names = FALSE)

# In excel generate a new list of genes that replaces the ENSEMBLID for the ENTREZID you obtained in lane 167.
# An example will be included in the data folder once the work is ready for publication.


# Loading data with ENTREZID + LFC 

df <- read.csv('./Data/list_of_genes_LFC_for_GO_ENTREZ.csv', header=TRUE)

head(df)

# name the vector
names(original_gene_list) <- df$ENTREZID.list

# omit any NA values (eliminate values not the names)
gene_list<-na.omit(original_gene_list)

# to omit the NA values and the refered gene name
gene_list <- gene_list[!is.na(names(gene_list))]

length(gene_list)

head(gene_list)



# GO over-representation analysis 
## Function --------------------------------------------------------------------


GO_results <- enrichGO(gene          = names(gene_list),
                       OrgDb         = "org.Hs.eg.db",
                       keyType       = "ENTREZID",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       minGSSize = 20, 
                       maxGSSize = 200, 
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

head(GO_results)

write.csv(GO_results,'./Results/GO_result.csv')

## Output ----------------------------------------------------------------------

#Dotplot:
#dotplot will show you the GO categories classified in your list of interest

dotp <- dotplot(GO_results, 
                color = "p.adjust",
                showCategory = 10,
                size = NULL,
                split = NULL,
                font.size = 10,
                orderBy = "x",
                label_format = 30,)


save_dotplot_pdf <- function(x, filename, width=10, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Save the file in pdf
save_dotplot_pdf(dotp, './Results/GO_result.pdf')




