####################################################################################
# TITLE - RNA-seq. Synergy analysis                                                #
# AUTHOR - Vincenzo Lagani and Mukhtar Sadykov                                     #
# DATE - December 20th, 2022                                                       #
# DESCRIPTION - This R script analyzes the synergy of different H1 linker histones #
# based on the transcriptional output of single knockout (SKO) and double knockout #
# (DKO) conditions. The script takes as input the gene expression data from SKO and# 
# DKO samples and performs comparative analysis to evaluate the combined effect of #
# H1 linker histones on gene expression.                                           #
####################################################################################

############################## IMPORTANT NOTES BEFORE STARTING ###############################
# 1. A "Data" folder has been included into the repository, this will include the processed  #
#    and normalized gene expression data, such as count matrices and metadata, once the work #
#    is ready for publication.                                                               #
# 2. A "Results" folder has been included into the repository. We have used this path to have# 
#    all "write.csv and "save" commands for our different pipelines                          #
##############################################################################################

# Script for performing Synergy analysis

# set up
rm(list = ls())

experiment.title <- "H1s_raw"
pacman::p_load(limma, edgeR, pheatmap, RColorBrewer,
ggplot2, ggpubr, qvalue, plyr, wesanderson,
GSEABase, grid, scales, WebGestaltR, stringr,
tidyverse, ggalluvial, corrplot, reshape2)
source('functions.R')

# control panel
data_folder <- 'Data'
res_folder_prefix <- 'Results'
fdr_thresholds <- c(0.01, 0.05, 0.1, 0.15, 0.2)
cpm_threshold <- 0.17 # at least 10 counts for at least 100 samples
synergy_levels <- c("undefined_LD", "undefined_LP",
"more_down_inverted", "more_down",
"less_down", "same", "less_up",
"more_up", "more_up_inverted")
zissou <- as.character(wes_palette("Zissou1", # color palette for the synergy levels
length(synergy_levels) - 2,
type = "continuous"))
zissou <- c('grey80', 'grey50', zissou)
names(zissou) <- synergy_levels
synergy_levels_refined <- c("undefined_LD", "undefined_LD_DKO_not_DF",
"undefined_LP", "undefined_LP_DKO_not_DF",
"more_down_inverted", "more_down_inverted_DKO_not_DF",
"more_down", "more_down_DKO_not_DF",
"less_down",  "less_down_DKO_not_DF",
"same", "same_DKO_not_DF",
"same_discordant_SKOs", "same_DKO_not_DF_discordant_SKOs",
"less_up", "less_up_DKO_not_DF",
"more_up", "more_up_inverted",
"more_up_DKO_not_DF", "more_up_inverted_DKO_not_DF")
zissou_refined <- as.character(wes_palette("Zissou1", # color palette for the synergy levels
length(synergy_levels_refined) - 2,
type = "continuous"))
zissou_refined <- c('grey80', 'grey80',
'grey50', 'grey50',
zissou_refined)
names(zissou_refined) <- synergy_levels_refined


# list of combinatorial experiments
sko_list <- c('0', '1', '2', '3', '4', '5', 'X')
sko1_list = c("H1.0vWT",      "H1.0vWT",       "H1.0vWT",       "H1.0vWT",       "H1.0vWT",       "H1.1vWT",       "H1.1vWT",       "H1.2vWT",       "H1.3vWT",       "H1.3vWT",       "H1.4vWT",       "H1.XvWT",       "H1.XvWT",       "H1.XvWT",       "H1.XvWT",       "H1.XvWT")
sko2_list = c("H1.1vWT",      "H1.2vWT",       "H1.3vWT",       "H1.4vWT",       "H1.5vWT",       "H1.4vWT",       "H1.5vWT",       "H1.4vWT",       "H1.4vWT",       "H1.5vWT",       "H1.5vWT",       "H1.0vWT",       "H1.2vWT",       "H1.3vWT",       "H1.4vWT",       "H1.5vWT")
syn_list = c("SynSKO_0plus1", "SynSKO_0plus2", "SynSKO_0plus3", "SynSKO_0plus4", "SynSKO_0plus5", "SynSKO_1plus4", "SynSKO_1plus5", "SynSKO_2plus4", "SynSKO_3plus4", "SynSKO_3plus5", "SynSKO_4plus5", "SynSKO_Xplus0", "SynSKO_Xplus2", "SynSKO_Xplus3", "SynSKO_Xplus4", "SynSKO_Xplus5")
add_list = c("AddSKO_0plus1", "AddSKO_0plus2", "AddSKO_0plus3", "AddSKO_0plus4", "AddSKO_0plus5", "AddSKO_1plus4", "AddSKO_1plus5", "AddSKO_2plus4", "AddSKO_3plus4", "AddSKO_3plus5", "AddSKO_4plus5", "AddSKO_Xplus0", "AddSKO_Xplus2", "AddSKO_Xplus3", "AddSKO_Xplus4", "AddSKO_Xplus5")
com_list = c("dk0_1vsWT",     "dk0_2vsWT",     "dk0_3vsWT",     "dk0_4vsWT",     "dk0_5vsWT",     "dk1_4vsWT",     "dk1_5vsWT",     "dk2_4vsWT",     "dk3_4vsWT",     "dk3_5vsWT",     "dk4_5vsWT",     "dkX_0vsWT",     "dkX_2vsWT",     "dkX_3vsWT",     "dkX_4vsWT",     "dkX_5vsWT")
experiment.titles = c("SynSKO_0plus1", "SynSKO_0plus2", "SynSKO_0plus3", "SynSKO_0plus4", "SynSKO_0plus5", "SynSKO_1plus4", "SynSKO_1plus5", "SynSKO_2plus4", "SynSKO_3plus4", "SynSKO_3plus5", "SynSKO_4plus5", "SynSKO_Xplus0", "SynSKO_Xplus2", "SynSKO_Xplus3", "SynSKO_Xplus4", "SynSKO_Xplus5")

# loading data
counts <- read.csv(file.path(data_folder, "all_KOs.csv"), row.names = 1)
meta <- read.csv(file.path(data_folder, "meta_all.csv"), row.names = 1)
meta <- meta[match(colnames(counts), row.names(meta)),]
anno <- read.csv(file.path(data_folder, "anno.csv"))
rownames(anno) <- anno$ensembl
tmp <- read.csv(file.path(data_folder, "anno_gene_type.csv"))
names(tmp) <- c('ensembl', 'gene_type')
rownames(tmp) <- tmp$ensembl

# merging the two annotation files
anno <- merge(anno, tmp, all.x = TRUE, all.y = FALSE)
anno$gene_type[is.na(anno$gene_type)] <- 'other'
row.names(anno) <- anno$ensembl

# looping across fdr thresholds
fdr_threshold <- fdr_thresholds[3]
for(fdr_threshold in fdr_thresholds){
    
    # results folders
    res_folder <- paste0('results_', fdr_threshold)
    dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)
    qc_folder <- file.path(res_folder, 'QC')
    dir.create(qc_folder, showWarnings = FALSE, recursive = TRUE)
    limma_folder <- file.path(res_folder, 'limma')
    dir.create(limma_folder, showWarnings = FALSE, recursive = TRUE)
    
    # CPM vs count quality control plot
    pdf(file.path(qc_folder, paste0(experiment.title, "-1_cpm-counts.pdf")))
    plot(cpm(counts)[, 1], counts[, 1], ylim = c(0, 50), xlim = c(0, 3))
    abline(h = 10, col = "red")
    abline(v = 0.17, col = "red")
    dev.off()
    
    # discarding lowly expressed genes
    keep <- rowSums(cpm(counts[]) > 0.17) >= 100 # at least 10 counts in at least 100 samples
    gExpr <- counts[keep, ]
    
    # creating the DGList object y
    y <- DGEList(gExpr)
    y <- calcNormFactors(y)
    
    # adding the gene annotation to y
    anno <- anno[match(rownames(y), rownames(anno)), ]
    y$genes <- anno
    
    # MDS plot on the normalized expression data
    pdf(file.path(qc_folder, paste0(experiment.title, "-2_mds_bc.pdf")))
    for (i in 1:length(colnames(meta))){
        mds(y, meta[ ,i], colnames(meta)[i])
    }
    plotMDS(y)
    dev.off()
    
    # creating the model matrix
    design <- model.matrix(~ 0 + Tx, meta)
    colnames(design) <- gsub("Tx", "", colnames(design))
    
    # voom normalization
    v <- voom(y, design, plot = TRUE, save.plot = TRUE)
    
    # voom QC
    pdf(file.path(qc_folder, paste0(experiment.title, "-3_voom.pdf")))
    plot(v$voom.xy, type = "p", pch=20, cex=0.16,
    main = "voom: Mean-variance trend",
    xlab = "log2( count size + 0.5 )",
    ylab = "Sqrt( standard deviation )")
    lines(v$voom.line, col="red")
    dev.off()
    
    # inter-sample correlation based on Batch information
    correlations <- duplicateCorrelation(v, design, block = meta$Batch)
    correlations <- correlations$consensus.correlation
    
    # fitting the limma models with both fixed and random effects
    fit <- lmFit(v, design, correlation = correlations, block = meta$Batch)
    
    # creating the contrast matrix
    cont.matrix <- create_contrast_matrix(design)
    
    # plotting the contrast matrix to ensure all contrast are correct
    cont.p <- t(cont.matrix)
    h <- pheatmap(cont.p,
    display_numbers = T, number_format = "%.0f",
    breaks = seq(-3, 1, by = 0.5),
    color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(12),
    cluster_cols = F, cluster_rows = F)
    pdf(file.path(qc_folder, paste0(experiment.title, "-design.pdf")), height = 15, width = 9)
    print(h)
    dev.off()
    
    # fitting the contrast
    fit.cont <- contrasts.fit(fit, cont.matrix)
    fit.cont <- eBayes(fit.cont)
    
    # mean variance trend
    pdf(file.path(qc_folder, paste0(experiment.title, "-mean_variance_trend.pdf")), height = 15, width = 9)
    plotSA(fit.cont, main = "Final model: Mean-variance trend", ylab = "Sqrt( standard deviation )")
    dev.off()
    
    # multiple testing correction
    summa.fit = decideTests(fit.cont, adjust.method = "fdr")
    
    # extracting results from each contrast and writing them on file
    res.list <- list()
    for (i in 1:length(colnames(fit.cont$contrasts))){
        x <- topTable(fit.cont, coef = i, sort.by = "p", n = Inf, confint = T)
        res.list[[i]] <- x
        names(res.list)[i] <- colnames(fit.cont$contrasts)[i]
        write.csv(x, file.path(limma_folder, paste0(experiment.title, "_DEGs_",
        colnames(fit.cont$contrasts)[i], ".csv")))
    }
    
    # plotting the volcano plots for each contrast
    pdf(file.path(res_folder, paste0(experiment.title, "-4_volcano-md-plots.pdf")))
    par(mfrow = c(1, 2))
    for (i in 1:length(colnames(fit.cont$contrasts))){
        plotMD(fit.cont, coef = i, status = summa.fit[, i], values = c(-1, 1))
        volcanoplot(fit.cont, coef = i, highlight = 10,
        names = fit.cont$genes$Gene_name)
    }
    dev.off()
    par(mfrow = c(1, 1))
    
    # plotting the expression values for the 3 most DE genes for each contrast
    pdf(file.path(res_folder, paste0(experiment.title, "-5_top3-expression-plots.pdf")))
    for (i in 1:length(colnames(fit.cont$contrasts))){
        x <- topTable(fit.cont, coef = i, sort.by = "p", n = Inf)
        cat("  \n\n### Plotting",  colnames(fit.cont$contrasts)[i], "  \n\n")
        for (j in 1:3){
            deg <- as.character(x[j,"ensembl"])
            p <- qplot(meta$Tx, v$E[deg, ],
            geom = "boxplot", fill = meta$Tx,
            ylab = "Normalized expression", xlab = "group",
            main = paste0(j, ". DEG: ", as.character(x[j, "Gene_name"]))) +
            geom_jitter() +
            rotate_x_text(angle = 45) +
            theme_bw(14)+
            theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
            print(p)
        }
    }
    dev.off()
    
    # matrix standard error from all the contrasts
    SE <- sqrt(fit.cont$s2.post) * fit.cont$stdev.unscaled
    
    # average SE from all experimental contrast
    tmp <- colnames(SE)[c(grep('dk', colnames(SE), fixed = TRUE),
    grep('H1.', colnames(SE), fixed = TRUE))]
    meanSE <- mean(SE[, tmp])
    
    # preparing results from all experiments
    all_experiments_results <- data.frame(ensembl = v$genes$ensembl,
    matrix('', length(v$genes$ensembl),
    length(experiment.titles)),
    matrix('', length(v$genes$ensembl),
    length(experiment.titles)))
    colnames(all_experiments_results) <- c('ensembl', paste0(experiment.titles,
    '.magnitude.syn'),
    paste0(experiment.titles,
    '.FDR'))
    rownames(all_experiments_results) <- all_experiments_results$ensembl
    all_experiments_results_refined <- data.frame(ensembl = v$genes$ensembl,
    matrix('', length(v$genes$ensembl),
    length(experiment.titles)),
    matrix('', length(v$genes$ensembl),
    length(experiment.titles)))
    colnames(all_experiments_results_refined) <- c('ensembl', paste0(experiment.titles,
    '.magnitude.syn.refined'),
    paste0(experiment.titles,
    '.FDR'))
    rownames(all_experiments_results_refined) <- all_experiments_results_refined$ensembl
    
    
    # producing results for each synergy (combinatorial) experiment
    i <- 16
    for (i in 1:length(syn_list)) {
        
        # info
        print(syn_list[i])
        
        # creating results directory
        current_res_folder <- file.path(res_folder, experiment.titles[i])
        dir.create(current_res_folder, showWarnings = FALSE, recursive = TRUE)
        
        # merging single KOs, additive, combinatorial and synergy results
        cols <- c("ensembl", "logFC", "P.Value", "adj.P.Val")
        tmp1 <- merge(res.list[sko1_list[i]][[1]][, c("Gene_name", cols)],
        res.list[sko2_list[i]][[1]][, cols],
        by=c("ensembl"))
        colnames(tmp1) <- c("ensembl", "Gene_name",
        paste0(sko1_list[i], c('.logFC', '.pvalue', '.FDR')),
        paste0(sko2_list[i], c('.logFC', '.pvalue', '.FDR')))
        tmp2 <- merge(tmp1,
        res.list[add_list[i]][[1]][, cols],
        by=c("ensembl"))
        colnames(tmp2) <- c(colnames(tmp1),
        paste0('Additive', c('.logFC', '.pvalue', '.FDR')))
        tmp3 <- merge(tmp2,
        res.list[com_list[i]][[1]][, cols],
        by=c("ensembl"))
        colnames(tmp3) <- c(colnames(tmp2),
        paste0('Combinatorial', c('.logFC', '.pvalue', '.FDR')))
        log2FC.matrix <- merge(tmp3,
        res.list[syn_list[i]][[1]][, cols],
        by=c("ensembl"))
        colnames(log2FC.matrix) <- c(colnames(tmp3),
        paste0('Synergistic', c('.logFC', '.pvalue', '.FDR')))
        rownames(log2FC.matrix) <- log2FC.matrix$ensembl
        
        # adding the synergy categories
        log2FC.matrix <- categorize.synergy(log2FC.matrix, meanSE, fdr_threshold)
        log2FC.matrix$magnitude.syn <- factor(log2FC.matrix$magnitude.syn,
        levels = synergy_levels)
        
        # refining the synergy categories
        log2FC.matrix <- categorize.synergy.refined(log2FC.matrix, meanSE, fdr_threshold,
        sko1_list[i], sko2_list[i])
        log2FC.matrix$magnitude.syn.refined <- factor(log2FC.matrix$magnitude.syn.refined,
        levels = synergy_levels_refined)
        
        # writing the synergy table on file
        write.csv(log2FC.matrix,
        file.path(current_res_folder, paste0(experiment.title, "_LogFC-FDR-synergy_matrix_",
        experiment.titles[i], ".csv")))
        
        # adding to the all results matrix
        if(!all(rownames(log2FC.matrix) == rownames(all_experiments_results))){
            stop('Different row names')
        }
        all_experiments_results[[paste0(experiment.titles[i], '.magnitude.syn')]] <-
        log2FC.matrix$magnitude.syn
        all_experiments_results[[paste0(experiment.titles[i], '.FDR')]] <-
        log2FC.matrix$Combinatorial.FDR
        if(!all(rownames(log2FC.matrix) == rownames(all_experiments_results_refined))){
            stop('Different row names (refined)')
        }
        all_experiments_results_refined[[paste0(experiment.titles[i], '.magnitude.syn.refined')]] <-
        log2FC.matrix$magnitude.syn.refined
        all_experiments_results_refined[[paste0(experiment.titles[i], '.FDR')]] <-
        log2FC.matrix$Combinatorial.FDR
        
        # alluvial plot for the significance of additive (expected) effect
        toPlot <- data.frame(sko1_sign = ifelse(log2FC.matrix[[paste0(sko1_list[i], '.FDR')]] < fdr_threshold,
        'sign', 'not_sign'),
        add_sign = ifelse(log2FC.matrix[['Additive.FDR']] < fdr_threshold, 'sign', 'not_sign'),
        sko2_sign = ifelse(log2FC.matrix[[paste0(sko2_list[i], '.FDR')]] < fdr_threshold,
        'sign', 'not_sign'))
        toPlot <- toPlot %>% group_by(sko1_sign, add_sign, sko2_sign) %>% summarise(freq = n())
        
        # alluvial plot
        png(file.path(current_res_folder, paste0(experiment.title, "--alluvial_skos_additive_",
        experiment.titles[i], ".png")),
        width = 2800, height = 2800, res = 300)
        p <- ggplot(toPlot,
        aes(y = freq, axis1 = sko1_sign, axis2 = add_sign, axis3 = sko2_sign)) +
        geom_alluvium(aes(fill = add_sign), alpha = 0.7, reverse = FALSE) +
        geom_stratum(width = 1/2, reverse = FALSE) +
        geom_text(stat = "stratum",
        aes(label = after_stat(stratum)), reverse = FALSE) +
        scale_x_continuous(breaks = 1:3, labels = c("First SKO",
        "Estimated additive effect",
        "Second SKO")) +
        scale_fill_discrete(name = 'significance') +
        theme_bw()
        plot(p)
        dev.off()
        
        # alluvial plot for the significance of combinatorial (actual) effect
        toPlot <- data.frame(Additive = ifelse(log2FC.matrix[['Additive.FDR']] < fdr_threshold,
        'sign', 'not_sign'),
        Synergistic = ifelse(log2FC.matrix[['Synergistic.FDR']] < fdr_threshold, 'sign', 'not_sign'),
        Combinatorial = ifelse(log2FC.matrix[['Combinatorial.FDR']] < fdr_threshold,
        'sign', 'not_sign'))
        toPlot <- toPlot %>% group_by(Additive, Synergistic, Combinatorial) %>% summarise(freq = n())
        
        # alluvial plot
        png(file.path(current_res_folder, paste0(experiment.title, "--alluvial_synergistic_",
        experiment.titles[i], ".png")),
        width = 2800, height = 2800, res = 300)
        p <- ggplot(toPlot,
        aes(y = freq, axis1 = Additive, axis2 = Synergistic, axis3 = Combinatorial)) +
        geom_alluvium(aes(fill = Synergistic), alpha = 0.7, reverse = FALSE) +
        geom_stratum(width = 1/2, reverse = FALSE) +
        geom_text(stat = "stratum",
        aes(label = after_stat(stratum)), reverse = FALSE) +
        scale_x_continuous(breaks = 1:3, labels = c("Additive (estimated)",
        "Synergistic (difference)",
        "Combinatorial (actual)")) +
        scale_fill_discrete(name = 'significance') +
        theme_bw()
        plot(p)
        dev.off()
        
        # computing and writing statistics on the synergy categories
        tmp <- table(log2FC.matrix$magnitude.syn)
        genes.per.category <- data.frame(category = factor(names(tmp), levels = synergy_levels),
        freq = as.numeric(tmp))
        genes.per.category$percent <- paste0(round(genes.per.category$freq*100/sum(genes.per.category$freq), 0), " %")
        write.csv(genes.per.category, file.path(current_res_folder,
        paste0(experiment.title, "_gene-count_synergy-categories_",
        experiment.titles[i], ".csv")))
        
        # barplot
        q <- ggplot(genes.per.category, aes(x = category, y = freq, fill = category)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=zissou) +
        theme_bw() +
        theme(legend.position = 'none') +
        coord_flip()
        pdf(file.path(current_res_folder, paste0(experiment.title, "-8_synergy-categories_barplot_",
        experiment.titles[i], ".pdf")))
        print(q)
        dev.off()
        
        # computing and writing statistics on the synergy categories (refined)
        tmp <- table(log2FC.matrix$magnitude.syn.refined)
        genes.per.category <- data.frame(category = factor(names(tmp), levels = synergy_levels_refined),
        freq = as.numeric(tmp))
        genes.per.category$percent <- paste0(round(genes.per.category$freq*100/sum(genes.per.category$freq), 0), " %")
        write.csv(genes.per.category, file.path(current_res_folder,
        paste0(experiment.title, "_gene-count_synergy-categories_refined_",
        experiment.titles[i], ".csv")))
        
        # barplot (refined)
        q <- ggplot(genes.per.category, aes(x = category, y = freq, fill = category)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=zissou_refined) +
        theme_bw() +
        theme(legend.position = 'none') +
        coord_flip()
        pdf(file.path(current_res_folder, paste0(experiment.title, "-8_synergy-categories_refined_barplot_",
        experiment.titles[i], ".pdf")))
        print(q)
        dev.off()
        
        # selecting the data for the heatmap
        tmp <- log2FC.matrix[ , c(paste0(sko1_list[i], '.logFC'),
        paste0(sko2_list[i], '.logFC'),
        "Additive.logFC","Combinatorial.logFC")]
        tmp <- tmp[order(rownames(tmp)), ]
        row_ann <- log2FC.matrix[ , 'magnitude.syn', drop = FALSE]
        row_ann$gene_type <- anno[rownames(row_ann), 'gene_type']
        row_ann <- row_ann[order(rownames(row_ann)), ]
        if(!all(rownames(row_ann) == row.names(tmp))){
            stop('row annotation misaligned')
        }
        
        # marking only protein coding, lncRNA and pseudogenes
        row_ann$gene_type[grep('pseudogene', row_ann$gene_type)] <- 'pseudogene'
        row_ann$gene_type[!(row_ann$gene_type %in% c('protein_coding',
        'lncRNA', 'pseudogene'))] <- 'other'
        
        # excluding undefined
        toKeep <- !grepl('undefined', row_ann$magnitude.syn, fixed = TRUE)
        row_ann <- row_ann[toKeep, ]
        tmp <- tmp[toKeep, ]
        
        # ordering the rows
        new_order <- order(row_ann$magnitude.syn,
        #row_ann$gene_type,
        tmp$Combinatorial.logFC,
        decreasing = TRUE)
        tmp <- tmp[new_order, ]
        row_ann <- row_ann[new_order, , drop = FALSE]
        
        # determinging the number of breaks for the color of the heatmap
        max_abs_log2FC <- max(c(abs(tmp$Additive.logFC),
        abs(tmp$Combinatorial.logFC)))
        breaks <- c(-1 * max_abs_log2FC,
        -1*log(seq(from = 9, to = 1.01, by = -0.1), 20),
        log(seq(from = 1, to = 9, by = 0.1), 20),
        max_abs_log2FC)
        nbreaks <- length(breaks)
        color <- colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(nbreaks - 1)
        
        # specify colors for the row annotation
        ann_colors = list(
        magnitude.syn = zissou,
        gene_type = c(other = "#D95F02", pseudogene = "#7570B3",
        lncRNA = "#E7298A", protein_coding = "#66A61E")
        )
        
        pdf(file.path(current_res_folder, paste0(experiment.title,
        "-9_synergy-categories_heatmap_",
        experiment.titles[i], ".pdf")),
        width = 7, height = 14)
        h <- pheatmap(tmp,
        cluster_rows = FALSE,
        scale = 'none',
        #annotation_row = row_ann[], # XXX
        annotation_row = row_ann[ , 'magnitude.syn', drop = FALSE],
        annotation_colors = ann_colors,
        cellwidth = 30,
        border_color = NA,
        breaks=breaks,
        cluster_cols = F,
        show_rownames = F,
        color = color,
        main = "logFC expected vs. measured")
        dev.off()
        
        # selecting the data for the heatmap
        tmp <- log2FC.matrix[ , c(paste0(sko1_list[i], '.logFC'),
        paste0(sko2_list[i], '.logFC'),
        "Additive.logFC","Combinatorial.logFC")]
        tmp <- tmp[order(rownames(tmp)), ]
        row_ann <- log2FC.matrix[ , 'magnitude.syn.refined', drop = FALSE]
        row_ann$gene_type <- anno[rownames(row_ann), 'gene_type']
        row_ann <- row_ann[order(rownames(row_ann)), ]
        if(!all(rownames(row_ann) == row.names(tmp))){
            stop('row annotation misaligned')
        }
        
        # marking only protein coding, lncRNA and pseudogenes
        row_ann$gene_type[grep('pseudogene', row_ann$gene_type)] <- 'pseudogene'
        row_ann$gene_type[!(row_ann$gene_type %in% c('protein_coding',
        'lncRNA', 'pseudogene'))] <- 'other'
        
        # excluding undefined
        toKeep <- !grepl('undefined', row_ann$magnitude.syn.refined, fixed = TRUE)
        row_ann <- row_ann[toKeep, ]
        tmp <- tmp[toKeep, ]
        
        # ordering the rows
        new_order <- order(row_ann$magnitude.syn.refined,
        #row_ann$gene_type,
        tmp$Combinatorial.logFC,
        decreasing = TRUE)
        tmp <- tmp[new_order, ]
        row_ann <- row_ann[new_order, , drop = FALSE]
        
        # determinging the number of breaks for the color of the heatmap
        max_abs_log2FC <- max(c(abs(tmp$Additive.logFC),
        abs(tmp$Combinatorial.logFC)))
        breaks <- c(-1 * max_abs_log2FC,
        -1*log(seq(from = 9, to = 1.01, by = -0.1), 20),
        log(seq(from = 1, to = 9, by = 0.1), 20),
        max_abs_log2FC)
        nbreaks <- length(breaks)
        color <- colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(nbreaks - 1)
        
        # specify colors for the row annotation
        ann_colors = list(
        magnitude.syn = zissou_refined,
        gene_type = c(other = "#D95F02", pseudogene = "#7570B3",
        lncRNA = "#E7298A", protein_coding = "#66A61E")
        )
        
        pdf(file.path(current_res_folder, paste0(experiment.title,
        "-9_synergy-categories_refined_heatmap_",
        experiment.titles[i], ".pdf")),
        width = 7, height = 14)
        h <- pheatmap(tmp,
        cluster_rows = FALSE,
        scale = 'none',
        #annotation_row = row_ann[], # XXX
        annotation_row = row_ann[ , 'magnitude.syn.refined', drop = FALSE],
        annotation_colors = ann_colors,
        cellwidth = 30,
        border_color = NA,
        breaks=breaks,
        cluster_cols = F,
        show_rownames = F,
        color = color,
        main = "logFC expected vs. measured")
        dev.off()
        
    }
    
    # barplot across all experiments
    toPlot <- all_experiments_results[, grep('magnitude.syn',
    colnames(all_experiments_results),
    fixed = TRUE)]
    colnames(toPlot) <- gsub('SynSKO_', '', colnames(toPlot), fixed = TRUE)
    colnames(toPlot) <- gsub('plus', '_', colnames(toPlot), fixed = TRUE)
    colnames(toPlot) <- gsub('.magnitude.syn', '', colnames(toPlot), fixed = TRUE)
    toPlot <- as.data.frame(pivot_longer(toPlot, cols = everything(),
    names_to = 'experiment', values_to = 'type'))
    pdf(file.path(res_folder, paste0(experiment.title, "-all_experiments_barplot.pdf")),
    width = 9, height = 4)
    p <- ggplot(toPlot, aes(x = experiment, fill = type)) +
    geom_bar(width = 0.8) +
    scale_fill_manual(values = zissou) +
    theme_bw()
    plot(p)
    dev.off()
    
    # correlation plot
    toPlot <- dcast(toPlot, experiment ~ type)
    rownames(toPlot) <- toPlot$experiment
    rownames(toPlot) <- gsub('.refined', '', rownames(toPlot), fixed = 'TRUE')
    toPlot$experiment <- NULL
    pdf(file.path(res_folder, paste0(experiment.title, "-all_experiments_corrplot.pdf")),
    width = 9, height = 5)
    corrplot(t(as.matrix(toPlot)), method = "shade",
    #addCoef.col = 'grey30', number.digits = 2, number.cex = 0.5, addCoefasPercent = TRUE,
    tl.col = 'black', is.corr = FALSE,
    col = brewer.pal(n = 9, name = 'Blues'), cl.align.text = 'l')
    dev.off()
    
    # chi2 test for assessing influence of SKO on synergy type
    pvalue_matrix <- matrix(NA, length(sko_list), dim(toPlot)[2])
    rownames(pvalue_matrix) <- sko_list
    colnames(pvalue_matrix) <- colnames(toPlot)
    for(m in 1:length(sko_list)){
        for(n in 1:length(colnames(toPlot))){
            
            tmp1 <- toPlot[, n]
            tmp2 <- rowSums(toPlot[, -n])
            names(tmp1) <- names(tmp2)
            toTest <- matrix(NA, 2, 2)
            toTest[1,1] <- sum(tmp1[grepl(sko_list[m], names(tmp1), fixed = TRUE)])
            toTest[1,2] <- sum(tmp1[!grepl(sko_list[m], names(tmp1), fixed = TRUE)])
            toTest[2,1] <- sum(tmp2[grepl(sko_list[m], names(tmp2), fixed = TRUE)])
            toTest[2,2] <- sum(tmp2[!grepl(sko_list[m], names(tmp2), fixed = TRUE)])
            pvalue_matrix[m,n] <- chisq.test(toTest)$p.value
            
        }
    }
    
    pvalue_matrix <- pvalue_matrix * prod(dim(pvalue_matrix))
    pvalue_matrix[pvalue_matrix > 1] <- 1
    pvalue_matrix_sign <- ifelse(pvalue_matrix < fdr_threshold, 1, 0)
    
    pdf(file.path(res_folder, paste0(experiment.title, "-all_experiments_corrplot_sign.pdf")),
    width = 9, height = 5)
    corrplot(pvalue_matrix_sign, method = 'shade', is.corr = FALSE,
    col = c('darkblue', 'darkred'), addgrid.col = 'grey',
    tl.col = 'black', cl.pos = 'n')
    dev.off()
    
    # barplot across all experiments (refined)
    toPlot <- all_experiments_results_refined[, grep('magnitude.syn.refined',
    colnames(all_experiments_results_refined),
    fixed = TRUE)]
    colnames(toPlot) <- gsub('SynSKO_', '', colnames(toPlot), fixed = TRUE)
    colnames(toPlot) <- gsub('plus', '_', colnames(toPlot), fixed = TRUE)
    colnames(toPlot) <- gsub('.magnitude.syn', '', colnames(toPlot), fixed = TRUE)
    toPlot <- as.data.frame(pivot_longer(toPlot, cols = everything(),
    names_to = 'experiment', values_to = 'type'))
    toPlot$experiment <- gsub('.refined', '', toPlot$experiment, fixed = TRUE)
    pdf(file.path(res_folder, paste0(experiment.title, "-all_experiments_refined_barplot.pdf")),
    width = 12, height = 4)
    p <- ggplot(toPlot, aes(x = experiment, fill = type)) +
    geom_bar(width = 0.8) +
    scale_fill_manual(values = zissou_refined) +
    theme_bw()
    plot(p)
    dev.off()
    
    # correlation plot (refined)
    toPlot <- dcast(toPlot, experiment ~ type)
    rownames(toPlot) <- toPlot$experiment
    rownames(toPlot) <- gsub('.refined', '', rownames(toPlot), fixed = 'TRUE')
    toPlot$experiment <- NULL
    pdf(file.path(res_folder, paste0(experiment.title, "-all_experiments_refined_corrplot.pdf")),
    width = 9, height = 5)
    corrplot(t(as.matrix(toPlot)), method = "shade",
    #addCoef.col = 'grey30', number.digits = 2, number.cex = 0.5, addCoefasPercent = TRUE,
    tl.col = 'black', is.corr = FALSE,
    col = brewer.pal(n = 9, name = 'Blues'), cl.align.text = 'l')
    dev.off()
    
    # chi2 test for assessing influence of SKO on synergy type
    pvalue_matrix <- matrix(NA, length(sko_list), dim(toPlot)[2])
    rownames(pvalue_matrix) <- sko_list
    colnames(pvalue_matrix) <- colnames(toPlot)
    for(m in 1:length(sko_list)){
        for(n in 1:length(colnames(toPlot))){
            
            tmp1 <- toPlot[, n]
            tmp2 <- rowSums(toPlot[, -n])
            names(tmp1) <- names(tmp2)
            toTest <- matrix(NA, 2, 2)
            toTest[1,1] <- sum(tmp1[grepl(sko_list[m], names(tmp1), fixed = TRUE)])
            toTest[1,2] <- sum(tmp1[!grepl(sko_list[m], names(tmp1), fixed = TRUE)])
            toTest[2,1] <- sum(tmp2[grepl(sko_list[m], names(tmp2), fixed = TRUE)])
            toTest[2,2] <- sum(tmp2[!grepl(sko_list[m], names(tmp2), fixed = TRUE)])
            pvalue_matrix[m,n] <- chisq.test(toTest)$p.value
            
        }
    }
    
    pvalue_matrix <- pvalue_matrix * prod(dim(pvalue_matrix))
    pvalue_matrix[pvalue_matrix > 1] <- 1
    pvalue_matrix_sign <- ifelse(pvalue_matrix < fdr_threshold, 1, 0)
    
    pdf(file.path(res_folder, paste0(experiment.title, "-all_experiments_refined_corrplot_sign.pdf")),
    width = 9, height = 5)
    corrplot(pvalue_matrix_sign, method = 'shade', is.corr = FALSE,
    col = c('darkblue', 'darkred'), addgrid.col = 'grey',
    tl.col = 'black', cl.pos = 'n')
    dev.off()
    
}
