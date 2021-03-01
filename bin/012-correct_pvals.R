#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c("-i", "--de_results"),
                        type = "character",
                        default = "de_results.tsv.gz",
                        help = "Dataframe containing merged DE results. DF
                            needs to have the following columns: pvalue, gene"
  ),

  optparse::make_option(c("-g", "--grouping_cols"),
                        type = "character",
                        default = "de_method,formula_passed,coef_value",
                        help = "Comma-separated list of columns to use to
                        stratify the dataframe to correct pvals."
  ),

  optparse::make_option(c("-o", "--output_file"),
                        type = "character",
                        default = "de_results.tsv.gz",
                        help = "Output file."
  ),

  optparse::make_option(c("-v", "--verbose"),
                        action = "store_true",
                        default = TRUE,
                        help = ""
  )
)

parser <- optparse::OptionParser(
  usage = "%prog",
  option_list = optionList,
  description = paste0(
    "Corrects diffential results using BH."
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
################################################################################

################################ Functions #####################################

plot_celltype_hits <- function(df,
                               cell_type_column,
                               p_val_column,
                               fdr_threshold = 0.05) {
  df$celltype <- as.factor(df[[cell_type_column]])
  df$fdr <- df[[p_val_column]]
  df <- df %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(
      n_hits = sum(fdr < fdr_threshold)
    ) %>%
    dplyr::select(celltype, n_hits) %>%
    unique( . )
  
  plot <- ggplot2::ggplot(df, ggplot2::aes(x=celltype,
                                           y=n_hits,
                                           fill=celltype)) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::labs(x="Cell type",
                  y="Number of genes (FDR < 0.05)",
                  fill=NULL) +
    ggplot2::scale_y_continuous(trans="log10") + 
    ggplot2::theme(legend.position="none")
    return(plot)
}

################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file <- arguments$options$output_file

# Get arguments
if (verbose) {
  print("Reading in the data...")
}

de_results <- read.csv(arguments$options$de_results,
                       sep = "\t",
                       header = T)
groupings <- strsplit(arguments$options$grouping_cols,
                      split = ",",
                      fixed = T)[[1]]

if (verbose) {
  print(sprintf("Performing correction using columns %s as groupings...",
                paste(groupings, collapse = ", ")))
}

de_results <- de_results %>%
  dplyr::group_by_at(.vars = groupings) %>%
  dplyr::mutate(
    qvalue_bh_allcelltypes = p.adjust(pvalue, method="BH")
  )
de_results$bh_correction_grouping <- arguments$options$grouping_cols

## Save result
if (verbose) {
    print("Writing DE results...")
}
gz_file <- gzfile(output_file, "w", compression = 9)
write.table(x = de_results,
            file = gz_file,
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)
close(gz_file)

if (verbose) {
    print("Done.")
}

## Plot results
if (verbose) {
  print("Plotting hits per cell type...")
}

rez <- de_results %>%
  dplyr::group_by_at(.vars = groupings) %>%
  dplyr::group_split() %>%
  lapply(. , function(x) {
    grouping <- sapply(groupings, function(group) {
      paste0(group, "::", unique(x[[group]]))
    })
    grouping_str <- paste(grouping, collapse="__")
    
    plot <- plot_celltype_hits(x,
                               "cell_label",
                               "qvalue_bh_allcelltypes",
                               fdr_threshold = 0.05)
    
    ggsave(sprintf("%s-celltype_hits.png", grouping_str),
           plot = plot,
           device = "png",
           dpi = 320)
  })

if (verbose) {
  print("Done.")
}

################################################################################
