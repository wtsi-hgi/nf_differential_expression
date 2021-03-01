#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c("-i", "--input_file"),
                        type = "character",
                        default = "de_results.tsv.gz",
                        help = "File containing DE results. Needs to contain
                        the following columns:
                        pvalue, qvalue_bh_per_cell_type, log2fc, mean_counts,
                        mean_cp10k"
  ),

  optparse::make_option(c("-t", "--target_var"),
                        type = "character",
                        default = "coef_value",
                        help = "Column used to facet the tests."
  ),

  optparse::make_option(c("-s", "--sig_var"),
                        type = "character",
                        default = "qvalue_bh_percelltype",
                        help = "Column used to check for significance."
  ),

  optparse::make_option(c("-x", "--sig_threshold"),
                        type = "double",
                        default = 0.05,
                        help = "Threshold used to check for significance."
  ),

  optparse::make_option(c("-l", "--sig_label"),
                        type = "character",
                        default = "FDR <= 0.05",
                        help = "Threshold used to check for significance."
  ),

  optparse::make_option(c("-o", "--out_file"),
                        type = "character",
                        default = "",
                        help = "Base output name."
  ),

  optparse::make_option(c("-v", "--verbose"),
                        action = "store_true",
                        default = FALSE,
                        help = ""
  )
)

parser <- optparse::OptionParser(
  usage = "%prog",
  option_list = optionList,
  description = paste0(
    "Plots results from differential gene expression."
  )
)

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
  optionStrings <- character()
  for (item in parserObj@options) {
    optionStrings <- append(optionStrings,
                            c(item@short_flag, item@long_flag))
  }
  optionStrings
}

optStrings <- getOptionStrings(parser)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
################################################################################

######################## Required Packages #####################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
################################################################################

################################ Functions #####################################

plot_volcano_plot <- function(df,
                              fc_col,
                              p_val_col,
                              facet_var,
                              sig_label) {
  df$neg_log10 <- -log10(df[[p_val_col]])
  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x = fc_col,
                                                  y = "neg_log10",
                                                  color = "significant")) +
    ggplot2::geom_point(size = 0.5, alpha = 0.25) +
    ggplot2::labs(x = "log2(FC)",
                  y = "-log10(p-value)",
                  color = sig_label) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("#CF9400", "Black")) +
    ggplot2::facet_wrap(as.formula(paste("~", facet_var)), nrow = 1)
  return(plot)
}


plot_ma_plot <- function(df,
                         mean_expr_col,
                         fc_col,
                         facet_var,
                         sig_label) {

  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x=mean_expr_col,
                                                  y=fc_col,
                                                  color="significant")) +
    ggplot2::geom_point(size = 0.5, alpha = 0.25) +
    ggplot2::scale_x_continuous(
        trans = "log10",
        labels = scales::comma_format()
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(x="Mean Expression",
                  y="log2(FC)",
                  color=sig_label) +
    ggplot2::scale_color_manual(values=c("#CF9400", "Black")) +
    ggplot2::facet_wrap(as.formula(paste("~", facet_var)), ncol = 1)
  return(plot)
}

################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file_base <- arguments$options$out_file

# Read DGE results in
dge_results <- read.csv(arguments$options$input_file,
                        sep ="\t",
                        header=T)

if (verbose) {
  cat(sprintf("Plotting DGE results for %s...\n",
              arguments$options$input_file))
}

# First find significance
sig_var <- arguments$options$sig_var
threshold <- arguments$options$sig_threshold

dge_results <- dge_results[!is.na(dge_results[[sig_var]]),]
dge_results$significant <- apply(dge_results, 1, function(x) {
  return(abs(as.numeric(x[[sig_var]])) <= threshold)
})

dge_results$significant <- factor(dge_results$significant,
                                  levels = c(TRUE, FALSE),
                                  labels = c("True", "False"))

facet_var <- arguments$options$target_var
sig_label <- arguments$options$sig_label
## Plot results
vol_plot <- plot_volcano_plot(
    dge_results,
    'log2fc',
    'pvalue',
    facet_var,
    sig_label
)
# ggplot2::ggsave(sprintf("%s-plot_volcano.png", output_file_base),
#        plot = vol_plot,
#        device = "png",
#        dpi = 320)
# png(
#     file =sprintf("%s-plot_volcano.png", output_file_base),
#     height = 10, width = 12, units = "cm", res = 320
# )
#     print(vol_plot)
# dev.off()

ma_plot <- plot_ma_plot(
    dge_results,
    'mean_counts',
    'log2fc',
    facet_var,
    sig_label
)
# ggplot2::ggsave(sprintf("%s-plot_ma.png", output_file_base),
#        plot = ma_plot,
#        device = "png",
#        dpi = 320)
# png(
#    file =sprintf("%s-plot_ma.png", output_file_base),
#    height = 10, width = 12, units = "cm", res = 320
# )
#    print(ma_plot)
# dev.off()


if (verbose) {
  cat("Done.\n")
}

################################################################################
