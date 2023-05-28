# script containing ad-hoc functions used in the analysis

# Function for categorizing synergy categories
# res_matrix: limma results with the following columns:
# - Synergistic.logFC
# - Synergistic.FDR
# - Additive.logFC
# - Additive.FDR
# - Combinatorial.logFC
# - Combinatorial.FDR
# logFC_threshold: logFC threshold
# fdr_threshold: significance threshold
categorize.synergy <- function(res_matrix, logFC_threshold, fdr_threshold){

  # creating the additional column
  res_matrix$magnitude.syn <- ''

  # looping over genes
  for (i in 1:length(res_matrix$Gene_name)){

    # first classification: is the synergy logFC significant?
    # In other words: is the difference between the actual and the expected
    # effect any significant?
    if(res_matrix$Synergistic.FDR[i] < fdr_threshold){

      # in this case the synergy logFC is significant

      # is the synergy logFC also in the +/-logFC_threshold limits?
      if(abs(res_matrix$Synergistic.logFC[i]) >= logFC_threshold){

        # here we have that the synergistic effect is both significant and
        # beyond the +/- logFC_threshold
        # we define different categories

        # positive synergy cases
        if(res_matrix$Synergistic.logFC[i] >= logFC_threshold){

          # more_up_inverted
          if(res_matrix$Additive.logFC[i] < 0 &
             res_matrix$Combinatorial.logFC[i] >= 0){
            res_matrix$magnitude.syn[i] <- 'more_up_inverted'
          }

          # more_up
          if(res_matrix$Additive.logFC[i] >= 0){
            res_matrix$magnitude.syn[i] <- 'more_up'
          }

          # less_down
          if(res_matrix$Additive.logFC[i] < 0 &
             res_matrix$Combinatorial.logFC[i] < 0){
            res_matrix$magnitude.syn[i] <- 'less_down'
          }

        }

        # negative synergy cases
        if(res_matrix$Synergistic.logFC[i] < -1 * logFC_threshold){

          # more_down_inverted
          if(res_matrix$Additive.logFC[i] >= 0 &
             res_matrix$Combinatorial.logFC[i] < 0){
            res_matrix$magnitude.syn[i] <- 'more_down_inverted'
          }

          # more_down
          if(res_matrix$Additive.logFC[i] < 0){
            res_matrix$magnitude.syn[i] <- 'more_down'
          }

          # less_up
          if(res_matrix$Additive.logFC[i] >= 0 &
             res_matrix$Combinatorial.logFC[i] >= 0){
            res_matrix$magnitude.syn[i] <- 'less_up'
          }

        }


      }else{

        # this is a borderline case: synergy logFC is significantly different but
        # the differnce between actual and expected synergistic effect is not beyond
        # the logFC threshoold. We label it "undefined, low difference"
        res_matrix$magnitude.syn[i] <- 'undefined_LD'

      }

    }else{

      # in this case the synergy logFC is not significant

      # is the synergy logFC also in the +/-logFC_threshold limits?
      if(abs(res_matrix$Synergistic.logFC[i]) >= logFC_threshold){

        # this is a borderline case: synergy logFC is quite large but not significant
        # we label it as "undefined, low power"
        res_matrix$magnitude.syn[i] <- 'undefined_LP'

      }else{

        # here we have that expected and actual synergistic effect are the same,
        # within a +/- logFC_threshold interval
        res_matrix$magnitude.syn[i] <- 'same'

      }

    }

  }

  # return
  return(res_matrix)

}

# Function for refining synergy categories
# res_matrix: limma results with the following columns:
# - Synergistic.logFC
# - Synergistic.FDR
# - Additive.logFC
# - Additive.FDR
# - Combinatorial.logFC
# - Combinatorial.FDR
# - magnitude.syn
# logFC_threshold: logFC threshold
# fdr_threshold: significance threshold
# sko1: name of the first intervention
# sko2: name of the first intervention
categorize.synergy.refined <- function(res_matrix, logFC_threshold, fdr_threshold,
                                       sko1, sko2){

  # creating the additional column
  res_matrix$magnitude.syn.refined <- as.character(res_matrix$magnitude.syn)

  # looping over genes
  for (i in 1:length(res_matrix$Gene_name)){

    # first classification: is the DKO logFC significant?
    if(res_matrix$Combinatorial.FDR[i] >= fdr_threshold){
      res_matrix$magnitude.syn.refined[i] <-
        paste0(res_matrix$magnitude.syn.refined[i], '_DKO_not_DF')
    }

    # second classification, only for the "same"
    if(res_matrix$magnitude.syn[i] == 'same'){

      # SKO differentially expressed with discordant signs
      if(res_matrix[i, paste0(sko1, '.FDR')] < fdr_threshold &&
         res_matrix[i, paste0(sko2, '.FDR')] < fdr_threshold &&
         res_matrix[i, paste0(sko1, '.logFC')] *
         res_matrix[i, paste0(sko2, '.logFC')] < 0){
        res_matrix$magnitude.syn.refined[i] <-
          paste0(res_matrix$magnitude.syn.refined[i], '_discordant_SKOs')
      }

    }


  }

  # return
  return(res_matrix)

}

# matrix of contrasts
create_contrast_matrix <- function(design){

  cont.matrix <- makeContrasts(
    dk0_1vsWT = DKO0_1 - WT_dko,
    dk0_2vsWT = DKO0_2 - WT_dko,
    dk0_3vsWT = DKO0_3 - WT_dko,
    dk0_4vsWT = DKO0_4 - WT_dko,
    dk0_5vsWT = DKO0_5 - WT_dko,

    #dk1_2vsWT = DKO1_2 - WT_dko,
    dk1_4vsWT = DKO1_4 - WT_dko,
    dk1_5vsWT = DKO1_5 - WT_dko,

    #dk2_3vsWT = DKO2_3 - WT_dko,
    dk2_4vsWT = DKO2_4 - WT_dko,
   # dk2_5vsWT = DKO2_5 - WT_dko,

    #dk3_1vsWT = DKO3_1 - WT_dko,
    dk3_4vsWT = DKO3_4 - WT_dko,
    dk3_5vsWT = DKO3_5 - WT_dko,

    dk4_5vsWT = DKO4_5 - WT_dko,

    dkX_0vsWT = DKOX_0 - WT_dko,
    #dkX_1vsWT = DKOX_1 - WT_dko,
    dkX_2vsWT = DKOX_2 - WT_dko,
    dkX_3vsWT = DKOX_3 - WT_dko,
    dkX_4vsWT = DKOX_4 - WT_dko,
    dkX_5vsWT = DKOX_5 - WT_dko,

    H1.0vWT = SKO_0 - WT_sko,
    H1.1vWT = SKO_1 - WT_sko,
    H1.2vWT = SKO_2 - WT_sko,
    H1.3vWT = SKO_3 - WT_sko,
    H1.4vWT = SKO_4 - WT_sko,
    H1.5vWT = SKO_5 - WT_sko,
    H1.XvWT = SKO_X - WT_sko,

    AddSKO_0plus1 = SKO_0 + SKO_1 - 2*WT_sko, # expectation of depletion
    SynSKO_0plus1 = (DKO0_1 - WT_dko) - (SKO_0 - WT_sko) - (SKO_1 - WT_sko),

    AddSKO_0plus2 = SKO_0 + SKO_2 - 2*WT_sko, # expectation of depletion
    SynSKO_0plus2 = (DKO0_2 - WT_dko) - (SKO_0 - WT_sko) - (SKO_2 - WT_sko),

    AddSKO_0plus3 = SKO_0 + SKO_3 - 2*WT_sko, # expectation of depletion
    SynSKO_0plus3 = (DKO0_3 - WT_dko) - (SKO_0 - WT_sko) - (SKO_3 - WT_sko),

    AddSKO_0plus4 = SKO_0 + SKO_4 - 2*WT_sko, # expectation of depletion
    SynSKO_0plus4 = (DKO0_4 - WT_dko) - (SKO_0 - WT_sko) - (SKO_4 - WT_sko),

    AddSKO_0plus5 = SKO_0 + SKO_5 - 2*WT_sko, # expectation of depletion
    SynSKO_0plus5 = (DKO0_5 - WT_dko) - (SKO_0 - WT_sko) - (SKO_5 - WT_sko),

    #AddSKO_1plus2 = SKO_1 + SKO_2 - 2*WT_sko, # expectation of depletion
    #SynSKO_1plus2 = (DKO1_2 - WT_dko) - (SKO_1 - WT_sko) - (SKO_2 - WT_sko),

    AddSKO_1plus4 = SKO_1 + SKO_4 - 2*WT_sko, # expectation of depletion
    SynSKO_1plus4 = (DKO1_4 - WT_dko) - (SKO_1 - WT_sko) - (SKO_4 - WT_sko),

    AddSKO_1plus5 = SKO_1 + SKO_5 - 2*WT_sko, # expectation of depletion
    SynSKO_1plus5 = (DKO1_5 - WT_dko) - (SKO_1 - WT_sko) - (SKO_5 - WT_sko),

    #AddSKO_2plus3 = SKO_2 + SKO_3 - 2*WT_sko, # expectation of depletion
    #SynSKO_2plus3 = (DKO2_3 - WT_dko) - (SKO_2 - WT_sko) - (SKO_3 - WT_sko),

    AddSKO_2plus4 = SKO_2 + SKO_4 - 2*WT_sko, # expectation of depletion
    SynSKO_2plus4 = (DKO2_4 - WT_dko) - (SKO_2 - WT_sko) - (SKO_4 - WT_sko),

    #AddSKO_2plus5 = SKO_2 + SKO_5 - 2*WT_sko, # expectation of depletion
    #SynSKO_2plus5 = (DKO2_5 - WT_dko) - (SKO_2 - WT_sko) - (SKO_5 - WT_sko),

    #AddSKO_3plus1 = SKO_1 + SKO_3 - 2*WT_sko, # expectation of depletion
    #SynSKO_3plus1 = (DKO3_1 - WT_dko) - (SKO_3 - WT_sko) - (SKO_1 - WT_sko),

    AddSKO_3plus4 = SKO_4 + SKO_3 - 2*WT_sko, # expectation of depletion
    SynSKO_3plus4 = (DKO3_4 - WT_dko) - (SKO_3 - WT_sko) - (SKO_4 - WT_sko),

    AddSKO_3plus5 = SKO_5 + SKO_3 - 2*WT_sko, # expectation of depletion
    SynSKO_3plus5 = (DKO3_5 - WT_dko) - (SKO_3 - WT_sko) - (SKO_5 - WT_sko),

    AddSKO_4plus5 = SKO_5 + SKO_4 - 2*WT_sko, # expectation of depletion
    SynSKO_4plus5 = (DKO4_5 - WT_dko) - (SKO_4 - WT_sko) - (SKO_5 - WT_sko),

    AddSKO_Xplus0 = SKO_X + SKO_0 - 2*WT_sko, # expectation of depletion
    SynSKO_Xplus0 = (DKOX_0 - WT_dko) - (SKO_X - WT_sko) - (SKO_0 - WT_sko),

    #AddSKO_Xplus1 = SKO_X + SKO_1 - 2*WT_sko, # expectation of depletion
    #SynSKO_Xplus1 = (DKOX_1 - WT_dko) - (SKO_X - WT_sko) - (SKO_1 - WT_sko),

    AddSKO_Xplus2 = SKO_X + SKO_2 - 2*WT_sko, # expectation of depletion
    SynSKO_Xplus2 = (DKOX_2 - WT_dko) - (SKO_X - WT_sko) - (SKO_2 - WT_sko),

    AddSKO_Xplus3 = SKO_X + SKO_3 - 2*WT_sko, # expectation of depletion
    SynSKO_Xplus3 = (DKOX_3 - WT_dko) - (SKO_X - WT_sko) - (SKO_3 - WT_sko),

    AddSKO_Xplus4 = SKO_X + SKO_4 - 2*WT_sko, # expectation of depletion
    SynSKO_Xplus4 = (DKOX_4 - WT_dko) - (SKO_X - WT_sko) - (SKO_4 - WT_sko),

    AddSKO_Xplus5 = SKO_X + SKO_5 - 2*WT_sko, # expectation of depletion
    SynSKO_Xplus5 = (DKOX_5 - WT_dko) - (SKO_X - WT_sko) - (SKO_5 - WT_sko),

    levels = design)

  return(cont.matrix )

}


# Function for categorizing synergy categories (original version)
# logFCmatrix: limma results with the following columns:
# - Gene_name
# - Synergistic.logFC
# - Synergistic.FDR
# - Additive.logFC
# - Additive.FDR
# - Combinatorial.logFC
# - Combinatorial.FDR
# meanSE: logFC threshold
#
categorize.synergy.old <- function(logFCmatrix, meanSE){

  # copying the original matrix and initializing the synergy magnitude
  m <- logFCmatrix
  m$magnitude.syn <- NA

  # looping over all genes
  for (i in 1:length(m$Gene_name)){

    # the synergy is above the meanSE threshold
    if (m$Synergistic.logFC[i] > meanSE){

      # the expected additive effect is below -1*threshold
      if (m$Additive.logFC[i] < -meanSE){

        # the actual combinatorial effect is above the threshold
        if (m$Combinatorial.logFC[i] > meanSE){
          # the expected is below the -1*threshold, the actual above threshold
          m$magnitude.syn[i] = "more.up"
        } else {
          # the expected is below the -1*threshold, the actual also below the -1*threshold,
          # but since actual - expected > meanSE, then less down
          m$magnitude.syn[i] = "less.down"
        }

      } else {
        # the expected additive effect is above the threshold

        # the sinergistic is positive, the expected additive is positive,
        # so the actual combinatorial must be more up
        m$magnitude.syn[i] = "more.up"
      }

    }else{

      # is the synergy below -meanSE threshold?
      if (m$Synergistic.logFC[i] < -meanSE){
        if (m$Additive.logFC[i] > meanSE){
          if (m$Combinatorial.logFC[i] < -meanSE){
            m$magnitude.syn[i] = "more.down"
          } else m$magnitude.syn[i] = "less.up"
        } else m$magnitude.syn[i] = "more.down"
      } else m$magnitude.syn[i] = "same"

    }

  }# end for

  # making the categry a factor and returning
  m$magnitude.syn <- as.factor(m$magnitude.syn)
  return(m)

}



# multi dimensional scaling plot on normalized count data
# normDGE: a normalized DGEList object
# metacol: a column from the meta data for coloring the points
mds <- function(normDGE, metacol, title){

  # specifying colors
  mcol <- as.factor(metacol)
  col <- rainbow(length(levels(mcol)), 1, 0.8, alpha = 0.5)[mcol]

  # plot, legend and title
  plotMDS(normDGE, col = col, pch = 16, cex = 2)
  legend("center",
         fill = rainbow(length(levels(mcol)), 1, 0.8),
         legend = levels(mcol),
         horiz = F,
         bty = "o",
         box.col="grey",
         xpd=TRUE)
  title(main=title)

}


cameraplusplots <- function(contrast, genesetlist, vobject, design, catcolors, title){
  tmp.list <- list()
  cam <- data.frame(matrix(ncol = 5, nrow = 0))
  for (i in 1:length(genesetlist)){
    cam.s <- camera(vobject, genesetlist[[i]], design, contrast = contrast, inter.gene.cor = 0.01)
    tmp.list[[i]] <- cam.s
    names(tmp.list)[i] <- names(genesetlist)[i]
    tmp.list[[i]]$category <- names(tmp.list[i])
    colnames(cam) <- names(tmp.list[[1]])
    cam <- rbind.data.frame(cam, tmp.list[[i]])
    print(paste0("Gene set categories run: ", i))
  }
  cam$neglogFDR <- -log10(cam$FDR)
  ## for plotting purposes only:
  cam$dirNeglogFDR <- cam$neglogFDR
  cam[(cam$Direction == "Down"), "dirNeglogFDR"] <- -cam[(cam$Direction == "Down"), "neglogFDR"]
  grob <- grobTree(textGrob(c("UP","DOWN"), x = c(0.94, 0.89), y = c(0.95, 0.05), hjust = 0, gp = gpar(fontsize = 13)))
  q <- ggplot(aes(x = cam$category, y = dirNeglogFDR, color = category), data = cam) +
    scale_color_manual(values = catcolors) +
    geom_jitter(aes(size = NGenes, alpha = neglogFDR), pch = 19, show.legend = F) +
    scale_size_continuous(range = c(4,16)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    geom_hline(yintercept = c(-1.3, 1.3), color = "red", alpha = 0.5) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-10, 10), oob = squish, labels = abs) +
    labs(x = "Gene set categories", y = "-log10(FDR)", title = title) +
    theme_bw(14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    annotation_custom(grob)
  print(q)
  cam$geneSet <- row.names(cam)
  cam10 <- as.data.frame(cam)
  cam10 <- cam10[order(cam10$FDR),]
  cam10 <- cam10[1:10,]
  grob <- grobTree(textGrob(c("DOWN","UP"), x = c(0.03, 0.9),  y=c(0.025), hjust = 0, gp = gpar(fontsize = 9, col = "grey60")))
  g <- ggplot(aes(x = geneSet, y = dirNeglogFDR, fill = category), data = cam10) +
    geom_col()+
    aes(reorder(stringr::str_wrap(geneSet, 60),-FDR), dirNeglogFDR) +
    xlab(NULL) +
    geom_hline(yintercept = c(-1.3, 1.3), color = "red", alpha = 0.3) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-10, 10), oob = squish, labels = abs) +
    labs(y = "-log10(FDR)", title = title) +
    scale_fill_manual(values = catcolors) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) +
    annotation_custom(grob)
  print(g)
  return(cam)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
power.compare.logFC <- function( sig1, sig2, N, N_other = c(2,4,6,8,10), alpha = 0.05, n_tests = 20000){
  d <- seq(0, 3, length.out=1000)
  alpha_multiple <- alpha / n_tests
  df <- lapply( N_other/N, function(n_scale){
    sigSq <- (sig1^2 + sig2^2) / n_scale
    cutoff <- qnorm( alpha_multiple/2, 0, sd = sqrt(sigSq), lower.tail = FALSE)
    p1 <- pnorm(-1*cutoff, d, sqrt(sigSq))
    p2 <- 1-pnorm(cutoff, d, sqrt(sigSq))
    data.frame(n_scale, d, power=p1+p2)
  })
  df <- do.call("rbind", df)
  ggplot(df, aes(d, power, color = as.factor(n_scale*N))) +
    geom_line() +
    theme_bw(14) +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    ylim(0, 1) +
    scale_color_discrete("Samples") +
    xlab(bquote(abs(logFC[observed] - logFC[expected]))) +
    ggtitle("Power versus difference in logFC")
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
stratify.by.syn.cat <- function(log2FC.matrix.sub){
  synergy.cat.list <- list(
    "less.down" = as.character(log2FC.matrix.sub[
    log2FC.matrix.sub$magnitude.syn == "less.down", "Gene_name"]),
    "less.up" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "less.up", "Gene_name"]),
    "more.down" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "more.down", "Gene_name"]),
    "more.up" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "more.up", "Gene_name"]),
    "same" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "same", "Gene_name"])
  )
  return(synergy.cat.list)
}


