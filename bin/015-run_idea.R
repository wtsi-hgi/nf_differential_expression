#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c('-i', '--de_results'),
                        type = 'character',
                        default = 'de_file',
                        help = 'TSV containing DE summary statistics. The
                        following columns are required: log2fc, pvalue, 
                        gene_symbol'
  ),
  
  optparse::make_option(c('-g', '--group_var'),
                        type = 'character',
                        default = 'coef_value',
                        help = 'Variable to group DE results by.'
  ),
  
  optparse::make_option(c('-b', '--max_beta_var'),
                        type = 'double',
                        default = 100,
                        help = 'Maximum effect size variance.'
  ),
  
  optparse::make_option(c('-c', '--min_cr'),
                        type = 'double',
                        default = 0.0025,
                        help = 'Minimum coverate rate for gene sets'
  ),
  
  optparse::make_option(c('-l', '--beta'),
                        type = 'integer',
                        default = NULL,
                        help = 'Initial value for the effect size in the 
                        MCMC algorithm.'
  ),
  
  optparse::make_option(c('-t', '--tau'),
                        type = 'character',
                        default = '(-2,0.5)',
                        help = 'Initial value for the coefficient of 
                        annotations/gene sets for the EM procedure.'
  ),
  
  optparse::make_option(c('-e', '--em_iter'),
                        type = 'integer',
                        default = 15,
                        help = 'Maximum iteration of EM algorithm.'
  ),
  
  optparse::make_option(c('-m', '--mcmc_iter'),
                        type = 'integer',
                        default = 1000,
                        help = 'Maximum iteration of MCMC algorithm.'
  ),
  
  optparse::make_option(c('-f', '--fit_toler'),
                        type = 'double',
                        default = 1e-5,
                        help = 'Tolerance for fitting the model.'
  ),
  
  optparse::make_option(c('-z', '--fdr_n_permutations'),
                        type = 'integer',
                        default = 10,
                        help = 'Number of permutations to run for FDR.'
  ),
  
  optparse::make_option(c('-d', '--database'),
                        type = 'character',
                        default = 'all',
                        help = 'Comma-separated list of databases to use from
                        the iDEA humanGeneSets compilation. Descriptions for 
                        databas can be found: 
                        https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
                        
                        Options:
                        Key Name        Description
                        c2.cgp          Chemical and genetic perturbations
                        c2.cp.biocarta  BioCarta
                        c2.cp.kegg      KEGG
                        c2.cp.reactome  Reactome
                        c2.cp           PID
                        c5.bp           GO biological process
                        c5.cc           GO cellular component
                        c5.mf           GO molecular function
                        c6.all          Oncogenic signatures
                        c7.all          Immunologic signatures
                        all             All gene sets'
  ),
  
  optparse::make_option(c('-o', '--output_file'),
                        type = 'character',
                        default = '',
                        help = 'Base output name.'
  ),
  
  optparse::make_option(c('-n', '--n_cores'),
                        type = 'integer',
                        default = 1,
                        help = 'Number of cores to use.'
  ),
  
  optparse::make_option(c('-v', '--verbose'),
                        action = 'store_true',
                        default = TRUE,
                        help = ''
  )
)

parser <- optparse::OptionParser(
  usage = '%prog',
  option_list = optionList,
  description = paste0(
    'Calculates differentially expressed genes using MAST.'
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
suppressPackageStartupMessages(library(iDEA))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggrepel))
################################################################################

################################ Functions #####################################

calculate_se <- function(df) {
  ## Get standard error for coefficient. Since DE methods don't output this,
  ## calculating using the formula found in 
  ## https://www.nature.com/articles/s41467-020-15298-6#Sec9
  return(abs(df$log2fc / qnorm(df$pvalue / 2.0, lower.tail = F)))
}

get_unique_symbols <- function(df) {
  return(df %>%
           group_by(gene_symbol) %>%
           group_split() %>%
           lapply( . , function(x) {
             ## If more than one unique value, get random and return
             if (nrow(x) > 1) {
               print(sprintf(paste('Gene symbol %s has more than 1 Ensembl ID.',
                                   'Selecting one at random to carry through.'),
                             x$gene_symbol[1]))
               index <- sample(1:nrow(x), 1)
               return(x[index,])
             }
             return(x)
             }) %>%
           do.call(rbind, .)
  )
}

plot_volcano_plot <- function(df,
                              coef_col,
                              p_val_col,
                              size_var,
                              cat_var,
                              facet_var,
                              coef_threshold = .8,
                              p_val_threshold = 0.05,
                              label_col) {
  df <- df[!(is.na(df[[coef_col]]))  & !(is.na(df[[p_val_col]])),]
  df$significant <- apply(df, 1, function(x) {
    return(abs(as.numeric(x[[coef_col]])) >= coef_threshold &
             abs(as.numeric(x[[p_val_col]])) <= p_val_threshold)
  })
  df$significant <- factor(df$significant,
                           levels = c(TRUE, FALSE),
                           labels = c('True', 'False'))
  
  df$neg_log10 <- -log10(df[[p_val_col]])
  plot <- ggplot2::ggplot(df,
                          ggplot2::aes_string(x = coef_col,
                                              y = 'neg_log10',
                                              color = 'significant')
                          ) +
    ggplot2::geom_point(shape = 19,
                        ggplot2::aes_string(fill = 'significant', 
                                            size = size_var,
                                            shape = 19),
                        alpha=0.7) +
    ggplot2::scale_radius() +
    ggplot2::labs(x = 'Enrichment',
                  y = expression(
                    paste(bold(-log[10]),
                          bold('('),
                          bolditalic(p),
                          bold('-value)')
                    )
                  ),
                  color = sprintf('Enrichment>=%s and FDR<=%s',
                                  round(coef_threshold, digits = 2),
                                  p_val_threshold),
                  size = 'Gene set size'
                  ) +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 1, 1, 'cm'),
                   plot.title = ggplot2::element_text(lineheight=.8, 
                                                      face='bold'),
                   axis.text = ggplot2::element_text(size = 24),
                   axis.line = ggplot2::element_line(colour = 'black'),
                   axis.ticks = ggplot2::element_line(colour = 'grey80'),
                   axis.title = ggplot2::element_text(size = 30,
                                                      face = 'bold'),
                   axis.title.y = ggplot2::element_text(margin = margin(t = 0,
                                                                        r = 20,
                                                                        b = 0,
                                                                        l = 0)
                   ),
                   legend.title = ggplot2::element_text(size=24,
                                                        face = 'bold'),
                   legend.text = ggplot2::element_text(size=24),
                   panel.background = ggplot2::element_rect(fill = NA,
                                                            color = "black"),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = 'gray'),
                   plot.background = ggplot2::element_blank(),
                   legend.key = ggplot2::element_rect(color = 'transparent',
                                                      fill = 'transparent'),
                   strip.text = ggplot2::element_text(size = 20)
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_val_threshold),
                        col = 'salmon',
                        linetype = 1,
                        size=1) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1,
                                                 override.aes = list(alpha = 1,
                                                                     shape = 21)
                                                 ),
                    color = ggplot2::guide_legend(order = 2,
                                                  override.aes=list(alpha = 1,
                                                                    shape = 21)
                                                  ),
                    fill = FALSE) +
    ggplot2::scale_color_manual(values = c('True'='salmon', 
                                           'False'='black')) +
    ggplot2::theme(legend.direction = 'horizontal',
                   legend.position = 'bottom',
                   legend.box = 'vertical',
                   legend.title.align = 0) +
    ggplot2::facet_grid(as.formula(paste(cat_var, '~', facet_var)))
  
  ## Get most significant terms and plot if there are any
  sig_terms <- df[order(df$neg_log10, decreasing = T),]
  sig_terms <- sig_terms[which(sig_terms$significant == "True"),]
  if (nrow(sig_terms) > 0) {
    sig_terms <- sig_terms[1:min(10, nrow(sig_terms)),]
    plot <- plot + 
      ggrepel::geom_text_repel(data = sig_terms,
                               ggplot2::aes_string(label = label_col),
                               force=1.0,
                               point.padding = ggplot2::unit(1.1,'lines'),
                               box.padding = ggplot2::unit(0.1, 'lines'),
                               size = 5,
                               col = 'black')
  }
  
  return(plot)
}

plot_boxplot <- function(df,
                         coef_col,
                         annotation_col,
                         sig_col,
                         sig_label,
                         facet_var,
                         sig_threshold = 0.05) {
  annotation_retain <- unique(
    df[[annotation_col]][df[[sig_col]] <= sig_threshold]
  )
  df <- df[df[[annotation_col]] %in% annotation_retain,]
  
  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x = coef_col,
                                                  y = annotation_col,
                                                  fill = sig_col)) +
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::labs(x = 'Enrichment',
                  y = 'Gene set',
                  fill = sig_label) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(as.formula(paste('~', facet_var)), nrow = 1) +
    ggplot2::theme(text = ggplot2::element_text(size=7))
  return(plot)
}

plot_bubble_plot <- function(df,
                             annot_id_col,
                             p_val,
                             color_col,
                             size_col,
                             facet_var) {
  cat_levels <- c('IMMUNOLOGIC SIGNATURES','CHEMICAL AND GENETIC PERTURBATIONS',
                 'GO BIOLOGICAL PROCESS','GO MOLECULAR FUNCTION',
                 'GO CELLULAR COMPONENT','ONCOGENIC SIGNATURES','REACTOME',
                 'KEGG','PID','BIOCARTA')
  df <- df[order(match(df[[color_col]], cat_levels)),]
  df$neg_log10 <- -log10(df[[p_val]])
  plot <- ggplot2::ggplot(df,
                          ggplot2::aes_string(x = annot_id_col,
                                              y = 'neg_log10',
                                              color = color_col)
                          ) +
    ggplot2::geom_point(shape = 19,
                        ggplot2::aes_string(fill = color_col, 
                                            size = size_col,
                                            shape = 19),
                        alpha=0.8) +
    ggplot2::scale_radius() +
    ggplot2::labs(x = 'Gene Sets',
                  y = expression(
                        paste(bold(-log[10]),
                              bold('('),
                              bolditalic(p),
                              bold('-value)')
                        )
                      )
                  ) +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 1, 1, 'cm'),
                   axis.text.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(lineheight=.8, 
                                                      face='bold'),
                   axis.text = ggplot2::element_text(size = 24),
                   axis.line = ggplot2::element_line(colour = 'black'),
                   axis.ticks = ggplot2::element_line(colour = 'grey80'),
                   axis.title = ggplot2::element_text(size = 30,
                                                      face = 'bold'),
                   axis.title.y = ggplot2::element_text(margin = margin(t = 0,
                                                                        r = 20,
                                                                        b = 0,
                                                                        l = 0)
                                                       ),
                   legend.title = ggplot2::element_text(size=24,
                                                        face = 'bold'),
                   legend.text = ggplot2::element_text(size=24),
                   panel.background = ggplot2::element_rect(fill = NA,
                                                            color = "black"),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = 'white'),
                   plot.background = ggplot2::element_blank(),
                   legend.key = ggplot2::element_rect(color = 'transparent',
                                                      fill = 'transparent'),
                   strip.text = ggplot2::element_text(size = 20)
                   ) +
    ggplot2::geom_hline(yintercept = 1.82,
                        col = 'black',
                        linetype = 2,
                        size=2) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1,
                                                 override.aes = list(alpha = 1,
                                                                     shape = 21)
                                                 ),
                    color = FALSE,
                    fill = FALSE) +
    ggplot2::labs(size = 'Gene set size',
                  color = 'Gene set category') +
    ggplot2::scale_color_manual(values = c('salmon','gold2','#42d4f4',
                                           '#3cb44b','chocolate2','#4363d8',
                                           '#bfef45','#911eb4','#f032e6',
                                           '#a9a9a9')) +
    ggplot2::scale_fill_manual(values = c('salmon','gold2','#42d4f4','#3cb44b',
                                          'chocolate2','#4363d8','#bfef45',
                                          '#911eb4','#f032e6','#a9a9a9')) +
    ggplot2::theme(legend.direction = 'horizontal',
                   legend.position = 'bottom',
                   legend.box = 'horizontal',
                   legend.title.align = 0) +
    ggplot2::facet_grid(as.formula(paste(facet_var, '~', color_col)),
                        scale = 'free_x',
                        space = 'free_x')
  #   ggrepel::geom_text_repel(data = Sig,
  #                            ggplot2::aes(label = Term),
  #                            force=1.0,
  #                            point.padding = ggplot2::unit(1.1,'lines'),
  #                            box.padding = ggplot2::unit(0.1, 'lines'),
  #                            vjust=-1.4,
  #                            hjust = 0.2,
  #                            size = 8,
  #                            direction='y',
  #                            nudge_x=0.2,
  #                            nudge_y = 0.2,
  #                            segment.size=0.2,
  #                            col = 'black')
  return(plot)
}

################################################################################


######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file_base <- arguments$options$output_file
n_cores <- arguments$options$n_cores
if (n_cores != 1) {
  n_cores <- n_cores - 1
}

# Read in data
if (verbose) {
  print('Reading in the data...')
}
de_genes <- read.csv(arguments$options$de_results,
                     sep='\t',
                     header=T)

# Get annotation data
# Currently using pre-compiled dataset from iDEA. See description:
# https://xzhoulab.github.io/iDEA/
# Subset according to `--database` args
data(humanGeneSets)
data(humanGeneSetsInfo)
rownames(humanGeneSetsInfo) <- humanGeneSetsInfo$gset ## Index by gset for later
if (arguments$options$database == 'all') {
  annot_data <- humanGeneSets
} else {
  if (verbose) {
    print(sprintf('Subsetting annotations to those apart of %s pathways...',
                  arguments$options$database))
  }
  databases <- strsplit(x = arguments$options$database,
                        split = ',',
                        fixed = T)[[1]]
  gs_ids <- humanGeneSetsInfo$gset[
    which(humanGeneSetsInfo$gsetName %in% databases)
  ]
  annot_data <- humanGeneSets[ , as.character(gs_ids)]
  if (verbose) {
    print(sprintf('After subsetting, there are %s gene sets remaining.',
                  ncol(annot_data)))
  }
}

# Get all arguments
max_beta_var <- arguments$options$max_beta_var
min_cover_rate <- arguments$options$min_cr
beta <- arguments$options$beta
tau <- as.double(
  strsplit(
    gsub(
      pattern='\\(|\\)',
      replacement='',
      x=arguments$options$tau,
      fixed=F
    ),
    split=',',
    fixed=T
  )[[1]]
)
em_iteration <- arguments$options$em_iter
mcmc_iteration <- arguments$options$mcmc_iter
fit_tolerance <- arguments$options$fit_toler
fdr_n <- arguments$options$fdr_n_permutations

## Create iDEA object and fit for each tested coefficient
group_vars <- strsplit(arguments$options$group_var,
                       split = ',',
                       fixed = T)[[1]]

iDEA_rl <- de_genes %>% 
  dplyr::group_by_at(.vars = group_vars) %>%
  dplyr::group_split() %>%
  lapply( . , function(x) {
    project <- paste(x[1, group_vars], collapse = ',')
    x$log2fc_se <- calculate_se(x)
    x <- get_unique_symbols(x)
    summary_data <- as.data.frame(x[ ,c('log2fc', 'log2fc_se')])
    rownames(summary_data) <- x$gene_symbol
    
    if (verbose) {
      print(sprintf('Calculating gene set enrichments for coefficient %s...',
                    project))
    }
    
    ## Build iDEA object
    if (verbose) {
      cat(sprintf(paste('Creating iDEA object using the arguments:\n',
                        'Maximum effect size variance: %s\n',
                        'Minimum coverage rate for gene sets: %s\n'),
                  max_beta_var,
                  min_cover_rate))
    }
    idea <- iDEA::CreateiDEAObject(summary_data,
                                   annot_data,
                                   project=project,
                                   max_var_beta=max_beta_var,
                                   min_precent_annot=min_cover_rate,
                                   num_core=n_cores)
    
    ## Fit model
    if (verbose) {
      beta_str <- beta
      if (is.null(beta)) {
        beta_str <- 'NULL'
      }
      cat(sprintf(paste('Fitting the iDEA model using the arguments:\n',
                        'Initial beta: %s\n',
                        'Initial tau: %s\n',
                        'EM iterations: %s\n',
                        'MCMC iterations: %s\n',
                        'Fit tolerance: %s\n'),
                  beta_str,
                  paste(tau, collapse = ', '),
                  em_iteration,
                  mcmc_iteration,
                  fit_tolerance))
    }
    idea <- iDEA::iDEA.fit(idea,
                           fit_noGS=FALSE,
                           init_beta=beta,
                           init_tau=tau,
                           min_degene=5,
                           em_iter=em_iteration,
                           mcmc_iter=mcmc_iteration,
                           fit.tol=fit_tolerance,
                           modelVariant = F,
                           verbose=verbose)
    
    ## Correct p-vals using louis method
    if (verbose) {
      print('Correcting the p-values...')
    }
    idea <- iDEA::iDEA.louis(idea)
    
    ## Also correct using BH FDR
    idea@gsea$qvalue_bh <- p.adjust(idea@gsea$pvalue, method = 'BH')
    
    # Now calculate FDR estimations using iDEAs permuted null distribution
    # method
    if (fdr_n != 0) {
      if (verbose) {
        print('Calculating calibrated FDR estimates...')
      }
      
      idea_null <- iDEA::iDEA.fit.null(idea,
                                       init_beta=beta,
                                       init_tau=tau,
                                       min_degene=5,
                                       em_iter=em_iteration,
                                       mcmc_iter=mcmc_iteration,
                                       fit.tol=fit_tolerance,
                                       modelVariant = F,
                                       numPermute = fdr_n,
                                       verbose=verbose)
      
      ## Correct p-vals using louis method
      idea.null <- iDEA::iDEA.louis(idea_null) 
      
      ## Calculate FDR
      df_FDR <- iDEA::iDEA.FDR(idea,
                               idea_null,
                               numPermute = fdr_n)
      rownames(df_FDR) <- df_FDR$annot_id
      idea@gsea$fdr_pvalue_permute <- df_FDR[idea@gsea$annot_id,]$fdr
    }
    
    ## Add shared run-specific information back to idea@gsea
    idea@gsea$gsea_method <- 'iDEA-zscore' ## in case we add beta effect model
    idea@gsea$de_method <- x$de_method[1]
    idea@gsea$reference_value <- x$reference_value[1]
    idea@gsea$coef_value <- idea@project # Add grouping var to gsea df 
    idea@gsea$formula_passed <- x$formula_passed[1]
    idea@gsea$formula <- x$formula[1]
    idea@gsea$cell_label_column <- x$cell_label_column[1]
    idea@gsea$cell_label <- x$cell_label[1]
    idea@gsea$max_beta_var <- max_beta_var
    idea@gsea$min_cover_rate <- min_cover_rate
    idea@gsea$beta_initial <- beta
    idea@gsea$tau_initial <- paste(tau, collapse = ',')
    idea@gsea$em_iter <- em_iteration
    idea@gsea$mcmc_iter <- mcmc_iteration
    idea@gsea$fit_tolerance <- fit_tolerance
    idea@gsea$fdr_n_permutations <- fdr_n
    
    ## Add categorical information for each id
    idea@gsea$gset_category_key <- humanGeneSetsInfo[idea@gsea$annot_id,
                                                     'gsetName']
    idea@gsea$gset_category <- humanGeneSetsInfo[idea@gsea$annot_id,
                                                 'gsetBioName']
    idea@gsea$gset_size <- Matrix::colSums(humanGeneSets[ , idea@gsea$annot_id])
    return(idea)
  })

gsea_rez <- do.call(rbind, lapply(iDEA_rl, function(x) {
  return(x@gsea) 
}))

# Plot results
vol_plot <- plot_volcano_plot(gsea_rez,
                              'annot_coef',
                              'pvalue_louis',
                              'gset_size',
                              'gset_category',
                              'coef_value',
                              0.8,
                              0.05,
                              'annot_id')
png(file = sprintf('%s-volcano_plot.png', output_file_base),
    width = 12500,
    height = 7500,
    res = 600
)
print(vol_plot)
dev.off()

# Plot results
bubble_plot <- plot_bubble_plot(gsea_rez,
                                'annot_id',
                                'pvalue_louis',
                                'gset_category',
                                'gset_size',
                                'coef_value')
png(file = sprintf('%s-bubble_plot.png', output_file_base),
    width = 12500,
    height = 7500,
    res = 600
)
print(bubble_plot)
dev.off()

if (fdr_n != 0) {
  ## Only do FDR -- pvals are inflated
  boxplot <- plot_boxplot(gsea_rez,
                          'annot_coef',
                          'annot_id',
                          'fdr_pvalue_permute',
                          'FDR',
                          'coef_value',
                          sig_threshold = 0.05)
  
  png(file = sprintf('%s-boxplot.png', output_file_base),
      units='px',
      width=1600,
      height=1600,
      res=300
  )
  print(boxplot)
  dev.off()
}

# Save mappings
if (verbose) {
  print('Writing GSEA results...')
}
gz_file <- gzfile(sprintf('%s-gsea_results.tsv.gz', output_file_base),
                  'w', 
                  compression = 9)
write.table(x=gsea_rez,
            file = gz_file,
            sep='\t',
            col.names=T,
            row.names=F,
            quote=F)
close(gz_file)

# Also save list of RDS idea files in case we ever need to access it
saveRDS(iDEA_rl,
        file=sprintf('%s-gsea_results.Rds.gz', output_file_base),
        compress=T)

if (verbose) {
  print('Done.')
}

################################################################################



