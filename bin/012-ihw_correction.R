#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c("-i", "--de_results"),
                        type = "character",
                        default = "de_results",
                        help = "Dataframe containing DE results. The dataframe
                            needs to have the following columns: pvalue, gene"
  ),

  optparse::make_option(c("-c", "--covariates"),
                        type = "character",
                        default = "cluster",
                        help = "Covariates to use in IHW correction. For
                            multiple covariates, needs to be a comma-separated
                            list."
  ),

  optparse::make_option(c("-a", "--alpha"),
                        type = "double",
                        default = 0.1,
                        help = "Alpha value to use. The alpha value is the
                            threshold under which one can reject the null
                            hypthesis."
  ),

  optparse::make_option(c("-g", "--grouping_cols"),
                        type = "character",
                        default = "de_method,formula_passed",
                        help = "Variables to group dataframes by before applying
                            correction method. Comma-separated list."
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
    "Corrects diffential results using IHW."
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
suppressPackageStartupMessages(library(IHW))
# suppressPackageStartupMessages(library(hexbin)) # for "decisionboundary" plot
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
################################################################################

################################ Functions #####################################

ihw_correct <- function(pval = "pvalue",
                        covariates,
                        df,
                        alpha_val = 0.1,
                        verbose = T) {
  if (length(covariates) == 1) {
    df$ihw_covariate <- df[[covariates]]
    df$ihw_covariate <- factor(df$ihw_covariate,
                               levels = df$ihw_covariate,
                               labels = df[[covariates]])
    if (verbose) {
      print(sprintf("Performing IHW correction using the formula: %s ~ %s",
                    pval,
                    covariates))
    }
  } else {
    if (verbose) {
      print(sprintf(paste("Calculating cartesion product of variables `%s`",
                          "to use as the IHW covariate."),
                    paste(covariates, collapse = ", ")))
    }
    df$ihw_covariate <- df %>%
      group_by_at(.vars = covariates) %>%
      group_indices()

    cov_labels <- apply(df[ , covariates ], 1, paste, collapse = "-" )
    df$ihw_covariate <- factor(df$ihw_covariate,
                               levels = df$ihw_covariate,
                               labels = cov_labels)
  }

  formula_str <- sprintf("%s ~ %s", pval, "ihw_covariate")
  formula <- formula(formula_str)
  rez <- IHW::ihw(formula,  data = df, alpha = alpha_val)
  return(rez)
}

plot_cov_qc <- function(
        df,
        pval = "pvalue",
        covariate,
        # max_bin = 20,
        grouping_str
    ) {
    # if (length(unique(df[[covariate]])) > max_bin) {
    #     if (verbose) {
    #         print(sprintf(
    #             "The covariate %s has more than %s levels, binning to %s.",
    #             covariate,
    #             max_bin,
    #             max_bin
    #         ))
    #     }
    if (is.numeric(df[[covariate]])) {
        df$cov <- factor(IHW::groups_by_filter(df[[covariate]], 20))
    } else {
        df$cov <- as.factor(df[[covariate]])
    }

    histo_plot <- ggplot2::ggplot(df, ggplot2::aes_string(x = pval))
    histo_plot <- histo_plot + ggplot2::geom_histogram(
        # binwidth = 0.025,
        boundary = 0
    )
    histo_plot <- histo_plot + ggplot2::facet_wrap(~ cov)
    histo_plot <- histo_plot + ggplot2::labs(x = "P-value")
    histo_plot <- histo_plot + ggplot2::theme_bw()
    histo_plot <- histo_plot + ggplot2::theme(
        axis.text.x = element_text(angle = 90)
    )
    png(
        file = sprintf("%s-%s-pvalue_histogram.png", grouping_str, covariate),
        units="px", width=1600, height=1600, res=300
    )
    print(histo_plot)
    dev.off()

    ecdf_plot <- ggplot2::ggplot(df, ggplot2::aes_string(
        x = pval,
        color = "cov"
    ))
    ecdf_plot <- ecdf_plot + ggplot2::stat_ecdf(geom = "step")
    ecdf_plot <- ecdf_plot + ggplot2::labs(x = "P-value", y = "ECDF")
    ecdf_plot <- ecdf_plot + ggplot2::theme_bw()
    if (length(unique(df$cov)) <= 8) {
        ecdf_plot <- ecdf_plot + ggplot2::scale_color_brewer(
            palette = "Dark2"
        )
    } else {
        ecdf_plot <- ecdf_plot + ggplot2::theme(legend.position = "none")
    }
    png(
        file = sprintf("%s-%s-pvalue_ecdf.png", grouping_str, covariate),
        units="px", width=1600, height=1600, res=300
    )
    print(ecdf_plot)
    dev.off()
    png(
        file = sprintf("%s-%s-pvalue_ecdf-facet.png", grouping_str, covariate),
        units="px", width=1600, height=1600, res=300
    )
    ecdf_plot <- ecdf_plot + ggplot2::facet_wrap(~ cov)
    ecdf_plot <- ecdf_plot + ggplot2::theme(
        axis.text.x = element_text(angle = 90)
    )
    print(ecdf_plot)
    dev.off()
    return(T)
}

################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file <- arguments$options$output_file

# Re-set options to allow multicore
# old <- options(stringsAsFactors = FALSE,
#                mc.cores=arguments$options$cores_available)
# on.exit(options(old), add = TRUE)

# Get arguments
if (verbose) {
  print("Reading in the data...")
}

de_results <- read.csv(
    arguments$options$de_results,
    sep = "\t",
    header = T
)
groupings <- strsplit(
    arguments$options$grouping_cols,
    split = ",",
    fixed = T
)[[1]]
covariates <- strsplit(
    arguments$options$covariates,
    split = ",",
    fixed = T
)[[1]]
alpha <- arguments$options$alpha

if (verbose) {
  print(sprintf("Performing correction using columns %s as groupings...",
                paste(groupings, collapse = ", ")))
}
corrected_list <- de_results %>%
  dplyr::group_by_at(.vars = groupings) %>%
  dplyr::group_split() %>%
  lapply(. , function(x) {
    grouping <- unique(x[,groupings])
    grouping_str <- paste(grouping, collapse="_")
    if (verbose) {
      print(sprintf("Correcting values for grouping: %s...",
                    paste(grouping, collapse = ", ")))
      print("Plotting QC metrics for each covariate...")
    }

    ## Generate QC plots
    for (cov in covariates) {
      if (verbose) {
        print(sprintf("Plotting QC metrics for covariate: %s...",
                      cov))
      }
      plot_cov_qc(x, "pvalue", cov, grouping_str)
    }

    ## Run IHW
    corrected_vals <- ihw_correct("pvalue", covariates, x, alpha, verbose)

    # Save raw results for use later
    if (verbose) {
      print("Writing IHW results...")
    }
    saveRDS(
        corrected_vals,
        file = sprintf("%s-ihw_results.Rds.gz", grouping_str),
        compress = T
    )
    corrected_vals_df <- as.data.frame(corrected_vals)
    corrected_vals_df$gene <- x$gene
    gz_file <- gzfile(
        sprintf("%s-ihw_results.tsv.gz", grouping_str),
        "w",
        compression = 9
    )
    write.table(
        x = corrected_vals_df,
        file = gz_file,
        sep = "\t",
        col.names = T,
        row.names = F,
        quote = F
    )
    close(gz_file)

    ## Save IHW plots
    pdf(
        file = sprintf("%s-ihw_results.pdf", grouping_str),
        height=5, width=6,
    )
        # TODO: In the case of a discrete covariate, change x axis to
        # the discrete covariate. If possible, sort x axis by the average
        # weight across folds.
        print(plot(corrected_vals))
        # NOTE: This plot does not make sense for a discrete covariate
        # print(plot(corrected_vals, what = "decisionboundary"))
    dev.off()
    png(
        file = sprintf("%s-ihw_results-pvalues.png", grouping_str),
        units = "px", width = 2600, height = 2600, res = 300
    )
    # Compare the expected vs the corrected p-values
    # NOTE: Corrected_vals_df$qvalue_ihw == corrected_vals_df$adj_pvalue
    corrected_vals_df$qvalue_ihw <- IHW::adj_pvalues(corrected_vals)
    plt <- ggplot2::ggplot(corrected_vals_df, ggplot2::aes(
        x = -log10(pvalue),
        y = -log10(adj_pvalue),
        color = "covariate"
    ))
    plt <- plt + ggplot2::theme_bw()
    plt <- plt + ggplot2::geom_point(alpha = 0.25)
    plt <- plt + ggplot2::geom_abline(slope = 1, intercept = 0)
    if (length(unique(corrected_vals_df$covariate)) > 8) {
        plt <- plt + ggplot2::facet_wrap(~ covariate, scales = "free")
        plt <- plt + ggplot2::theme(
            axis.text.x = element_text(angle = 90),
            legend.position = "none"
        )
    }
    print(plt)
    dev.off()

    ## Add back to df and return
    x$qvalue_ihw_allcelltypes <- IHW::adj_pvalues(corrected_vals)
    if (verbose) {
      print("Done.")
    }
    return(x)
  })

# Make the final dataframe
de_results <- do.call(rbind, corrected_list)
# Add covariates to the data.frame
de_results$ihw_covariates <- arguments$options$covariates

## Save result
if (verbose) {
    print("Writing DE results...")
}
gz_file <- gzfile(output_file, "w", compression = 9)
write.table(
    x = de_results,
    file = gz_file,
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
)
close(gz_file)

if (verbose) {
    print("Done.")
}

################################################################################
