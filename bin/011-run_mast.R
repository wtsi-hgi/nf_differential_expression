#!/usr/bin/env Rscript

######################## Required Packages #####################################
suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
################################################################################

################################ Functions #####################################

DE_get_input_type <- function() {
  return("logcp10")
}

DE_allow_random_effect <- function() {
  return(TRUE)
}

## Code is lifted from
## https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
DE_calculate_dge <- function(
  input_matrix,
  feature_metadata,
  sample_metadata,
  testing_var,
  coef_value,
  formula,
  method = "singlecell::bayesglm",
  n_cores = 1,
  verbose = TRUE
) {
  resolution <- strsplit(method, split="::", fixed=T)[[1]][1]
  if (resolution == "pseudobulk") {
    stop(paste("You are trying to run MAST with pseudobulk. MAST was designed",
               "for single-cell data"))
  }
  de_method <- strsplit(method, split="::", fixed=T)[[1]][2]
  
  if (any(grepl("|", attr(terms(formula), 'term.labels'), fixed=T))) {
    if (de_method != "glmer") {
      stop(sprintf(paste(
        "You are trying to run MAST with a random effect using method",
        "`%s`. MAST only runs with a random effect using glmer."
      )))
    }
  }

  if (n_cores > 1) {
    # Re-set options to allow multicore
    old <- options(stringsAsFactors = FALSE, mc.cores=n_cores)
    on.exit(options(old), add = TRUE)
  }

  # First, format data into sca object and fit model
  sca <- MAST::FromMatrix(
    exprsArray = input_matrix,
    cData = sample_metadata,
    fData = feature_metadata
  )

  # https://github.com/RGLab/MAST/issues/107
  ebayes = TRUE
  if (de_method == "glmer") {
    ebayes = FALSE
  }

  # MAST::zlm is the hurdle model
  # For more information: https://rdrr.io/bioc/MAST/man/zlm.html
  if (verbose) {
    cat("Fitting the model using MAST:zlm...\n")
  }
  zlm_fit <- MAST::zlm(
    formula,
    sca,
    method = de_method,
    ebayes = ebayes,
    silent = FALSE,
    parallel = TRUE
  )
  if (verbose) {
    cat("Done fitting the model.\n")
  }

  # Check for testing var in model
  if (!(coef_value %in% colnames(zlm_fit@coefD))) {
    cat(sprintf(paste("Returning empty dataframe for term %s",
                      "because it does't exist in the model. It was most",
                      "likely removed because it is linearly associated",
                      "with another term."),
                  coef_value))
    rez <- data.frame()
    return(rez)
  }

  if (verbose) {
    cat(sprintf("Performing Likelihood Ratio Test for the variable %s...\n",
                coef_value))
  }

  ## This is another way to calculate LRT
  ## Passing variable in as a character treats covariate as a continuous
  ## variable. Using summary(zlm_fit, doLRT=test_var) and adding logic
  ## to pass in variables correctly.
  # results <- MAST::lrTest(zlm_fit, test_var)
  results <- summary(zlm_fit, doLRT = coef_value, logFC = TRUE)
  results_dt <- results$datatable

  # Structure results into interpretable format
  fcHurdle <- merge(
    results_dt[
      contrast == coef_value & component =='H',
      .(primerid, `Pr(>Chisq)`)
    ], #hurdle P values
    results_dt[
      contrast == coef_value & component =='logFC',
      .(primerid, coef, ci.hi, ci.lo, z)
    ],
    by = 'primerid'
  ) #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle <- merge(fcHurdle, as.data.frame(mcols(sca)), by='primerid')

  ## Order by pvalue
  # fcHurdle <- fcHurdle[order(fcHurdle[[`Pr(>Chisq)`]], decreasing = F)]
  if (verbose) {
    cat("Done performing Likelihood Ratio Test.\n")
  }

  ## Rename MAST results for downstream analysis
  fcHurdle <- fcHurdle %>%
    dplyr::rename(
      gene = primerid,
      pvalue = `Pr(>Chisq)`,
      log2fc = coef,
      test_statistic = z ## Use test_statistic for downstream anaylsis
    ) %>%
    dplyr::arrange(pvalue) %>%
    as.data.frame()

  ## Add other necessary items
  fcHurdle$test_statistic_type <- "z_score"

  return(fcHurdle)
}

################################################################################
