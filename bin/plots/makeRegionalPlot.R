#!/usr/bin/env Rscript

makeRegionalPlotClass <- R6::R6Class(
  "makeRegionalPlotClass",
  public = list(
    plot.dat = NULL,
    recomb.dat = NULL,
    plot.name = NULL,
    ensembl = NULL,
    lead_var = NULL,
    snp_col = NULL,
    pval_col = NULL,
    initialize = function(plot.dat, recomb.dat, plot.name, ensembl, lead_var, snp_col, pval_col) {
      self$plot.dat <- plot.dat
      self$recomb.dat <- recomb.dat
      self$plot.name <- plot.name
      self$ensembl <- ensembl
      self$lead_var <- lead_var
      self$snp_col <- snp_col
      self$pval_col <- pval_col
    },
    plot = function() {
      plot.dat <- self$plot.dat
      recomb.dat <- self$recomb.dat
      plot.name <- self$plot.name
      ensembl <- self$ensembl
      lead_var <- self$lead_var
      snp_col <- self$snp_col
      pval_col <- self$pval_col

      plot.dat$pval <- subset(plot.dat, select = pval_col)[, 1]
      plot.dat$snpid <- subset(plot.dat, select = snp_col)[, 1]
      if (!is.null(lead_var)) {
        i_lead <- plot.dat$snpid == lead_var
      }

      ld_cat <- c(.2, .4, .6, .8, 1)
      col_cat <- c('darkblue', 'lightblue', 'green', 'orange', 'red')

      plot.dat$BP <- as.integer(plot.dat$BP)

      min.pos <- min(plot.dat$BP)
      max.pos <- max(plot.dat$BP)

      plot.dat$R2_cat <- 0
      for (i in 5:1) {
        plot.dat$R2_cat[plot.dat$R2 <= ld_cat[i]] <- i
      }

      chrom <- paste("chr", plot.dat$CHR[1], sep = "")
      chromstart <- min.pos
      chromend <- max.pos
      chrom_biomart <- gsub("chr", "", chrom)
      geneinfo <- getBM(attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "external_gene_name", "strand"),
                        filters = c("chromosome_name", "start", "end"),
                        values = list(chrom_biomart, chromstart, chromend),
                        mart = ensembl)
      names(geneinfo) <- c("chrom", "start", "stop", "gene", "strand")
      geneinfo$score <- "."
      geneinfo <- geneinfo[, c(1:4, 6, 5)]  
      r <- max.pos - min.pos
      geneinfo <- subset(geneinfo, start > (min.pos + r * .05) & start < (max.pos - r * .05))
      wdivider <- c(60, 50, 40) / 10
      s <- c(20, 40, 60)
      l <- length(unique(geneinfo$gene))
      wdivider <- wdivider[min(abs(l - s)) == abs(l - s)]
      wlayout <- ifelse(c(1:10) <= wdivider, 2, 1)

      keep.recomb <- subset(recomb.dat, recomb.dat[, 2] > min.pos & recomb.dat[, 3] < max.pos)

      png(paste0(plot.name, ".png"), height = 6, width = 6, units = 'in', res = 300)
      layout(matrix(wlayout, ncol = 1, byrow = TRUE))

      par(mar = c(0, 4.5, 2, 4.5))
      pg <- plotGenes(geneinfo, chrom, chromstart, chromend,
                      types = 'exon', bheight = .2,
                      plotgenetype = "box", bentline = F, packrow = T,
                      wigglefactor = 0.3, labelat = "middle", labeloffset = .5,
                      labeltext = T, col = "blue", xaxs = "i")

      par(mar = c(0, 4.5, 2, 4.5))
      plot(keep.recomb[, 2], keep.recomb[, 3], type = "l", col = "lightblue", lwd = 1, xlim = c(min.pos, max.pos),
           ylim = c(0, 100), xlab = "", ylab = "", axes = F, xaxs = "i")

      par(new = TRUE)
      plot(plot.dat$BP, -log10(plot.dat$pval),
           xlim = c(min.pos, max.pos),
           ylim = c(0, (-log10(plot.dat$pval[i_lead]) + 1)), type = "n", xlab = "", ylab = "", axes = F, xaxs = "i")

      for (i in 1:5) {
        ixx <- plot.dat$R2_cat == i
        points(plot.dat$BP[ixx], -log10(plot.dat$pval[ixx]), pch = 19, col = col_cat[i])
      }
      if (!is.null(lead_var)) {
        points(plot.dat$BP[i_lead], -log10(plot.dat$pval[i_lead]), pch = 17, col = "purple",
               cex = 1.5)
        text(plot.dat$BP[i_lead], (-log10(plot.dat$pval[i_lead]) + 1), labels = lead_var)
      }

      axis(2, las = 1)
      mtext(expression('-log'[10]~'(P)'), side = 2, at = ceiling(-log10(min(plot.dat$pval))) / 2, line = 2, cex = .75)

      ystep <- (-log10(plot.dat$pval[i_lead]) + 1) / 5

      axis(4, at = seq(0, ystep * 5, ystep), labels = seq(0, 100, 20), las = 1)
      mtext("Recombination rate (cM/Mb)", side = 4, at = ystep * 2.5, line = 2, cex = .75, las = 0)
      labelgenome(chrom, chromstart, chromend, n = 3, scale = "Mb", chromcex = .8, scalecex = .8, xaxs = "i")

      image.plot(legend.only = T, zlim = c(0, length(col_cat) - 1), col = col_cat, legend.shrink = .5, smallplot = c(.12, .14, .5, .8),
                 axis.args = list(at = seq(0.5, 3.5), labels = seq(.2, .8, .2)), legend.args = list(text = expression("r"^2), side = 3))

      dev.off()
    }
  )
)

# Wrapper for compatibility with existing pipeline calls
makeRegionalPlot <- function(plot.dat, recomb.dat, plot.name, ensembl, lead_var, snp_col, pval_col) {
  plotter <- makeRegionalPlotClass$new(plot.dat, recomb.dat, plot.name, ensembl, lead_var, snp_col, pval_col)
  plotter$plot()
}