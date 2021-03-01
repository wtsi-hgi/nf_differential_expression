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

  optparse::make_option(c('-r', '--ranking_var'),
                        type = 'character',
                        default = 'test_statistic',
                        help = 'Column of de_result to rank the genes by.'
  ),

  optparse::make_option(c('-q', '--sample_size'),
                        type = 'integer',
                        default = 101,
                        help = 'Sample size to use during error calculation.'
  ),

  optparse::make_option(c('-l', '--min_set_size'),
                        type = 'integer',
                        default = 1,
                        help = 'Minimum number of genes for a gene set.'
  ),

  optparse::make_option(c('-m', '--max_set_size'),
                        type = 'double',
                        default = Inf,
                        help = 'Maximum number of genes for a gene set.'
  ),

  optparse::make_option(c('-e', '--eps'),
                        type = 'double',
                        default = 1e-10,
                        help = 'Boundary for calculating the p-value.'
  ),

  optparse::make_option(c('-s', '--score_type'),
                        type = 'character',
                        default = 'std',
                        help = 'GSEA score type. Possible options: ("std",
                        "pos", "neg")'
  ),
  
  optparse::make_option(c('-t', '--unsigned_ranking'),
                        action = 'store_true',
                        default = FALSE,
                        help = 'If `ranking_var` should be signed.'
  ),

  optparse::make_option(c('-z', '--gsets_gene_matrix'),
                        type = 'character',
                        default = '',
                        help = '
                          TSV file containing gene sets with genes.
                        '
  ),
  
  optparse::make_option(c('-y', '--gsets_info_file'),
                        type = 'character',
                        default = '',
                        help = '
                          TSV file containing gene set info.
                        '
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
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggrepel))
################################################################################

################################ Functions #####################################

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
                  color = sprintf('Enrichment>=%s and P-value<=%s',
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

# Check if test statistic should be signed or not
if (arguments$options$unsigned_ranking) {
  print(paste(
    '`unsigned_ranking` flag is set. Taking absolute value of the ',
    '`ranking_var`.'
  ))
  de_genes[[arguments$options$ranking_var]] <- abs(
    de_genes[[arguments$options$ranking_var]]
  )
}

# Get annotation data
# Subset according to `--database` args
gene_set_genes <- read.csv(arguments$options$gsets_gene_matrix,
                           sep='\t',
                           header=T)
gene_set_info <- read.csv(arguments$options$gsets_info_file,
                           sep='\t',
                           header=T)
rownames(gene_set_info) <- gene_set_info$gene_set ## Index by gset for later
if (arguments$options$database == 'all') {
  annot_data <- gene_set_info
} else {
  if (verbose) {
    print(sprintf('Subsetting annotations to those apart of %s pathways...',
                  arguments$options$database))
  }
  databases <- strsplit(x = arguments$options$database,
                        split = ',',
                        fixed = T)[[1]]
  gs_ids <- gene_set_info$gene_set[
    which(gene_set_info$category %in% databases)
  ]
  annot_data <- gene_set_genes[ , as.character(gs_ids)]
  if (verbose) {
    print(sprintf('After subsetting, there are %s gene sets remaining.',
                  ncol(annot_data)))
  }
}

## Format annotation data to be compatible with fGSEA
annot_data <- sapply(colnames(annot_data), function(x) {
  return(list(
    rownames(annot_data)[which(annot_data[[x]] == 1)]
  ))
})

# Get all arguments
sample_size <- arguments$options$sample_size
min_set_size <- arguments$options$min_set_size
max_set_size <- arguments$options$max_set_size
eps <- arguments$options$eps
score_type <- arguments$options$score_type

## Run fGSEA for each tested coefficient
group_vars <- strsplit(arguments$options$group_var,
                       split = ',',
                       fixed = T)[[1]]

fgsea_rl <- de_genes %>%
  dplyr::group_by_at(.vars = group_vars) %>%
  dplyr::group_split() %>%
  lapply( . , function(x) {
    project <- paste(x[1, group_vars], collapse = ',')
    x <- get_unique_symbols(x)

    ## Rank genes by ranking_var
    summary_data <- x[[arguments$options$ranking_var]]
    names(summary_data) <- x$gene_symbol
    summary_data <- na.omit(summary_data)

    if (verbose) {
      print(sprintf('Calculating gene set enrichments for coefficient %s...',
                    project))
      cat(sprintf(paste('Performing fGSEA using the arguments:\n',
                        'Sample size parameter: %s\n',
                        'Maximum gene set size: %s\n',
                        'Minimum gene set size: %s\n',
                        'EPS value: %s\n',
                        'Score type: %s\n'),
                  sample_size,
                  max_set_size,
                  min_set_size,
                  eps,
                  score_type))
    }

    ## This will run fgsea.Multilevel, which gives more precise p-vals at a
    ## consequence for longer compute time. Could adjust to using fgsea.Simple
    ## See: https://github.com/ctlab/fgsea/blob/master/R/fgsea.R
    fgseaRes <- fgsea::fgseaMultilevel(pathways = annot_data,
                                       stats = summary_data,
                                       sampleSize = sample_size,
                                       eps = eps,
                                       minSize = min_set_size,
                                       maxSize = max_set_size,
                                       scoreType = score_type,
                                       nproc = n_cores,
                                       gseaParam = 1,
                                       BPPARAM = NULL)

    ## Rename fgsea results ot match iDEA
    fgseaRes <- fgseaRes %>%
      dplyr::rename(
        annot_id = pathway,
        annot_coef = NES,
        pvalue = pval,
        qvalue_bh = padj, ## fgsea automatically does BH correction
        n_gene_contained = size
      )

    ## Add shared run-specific information back to result
    fgseaRes$gsea_method <- 'fgsea-multi' ## in case we add beta effect model
    fgseaRes$de_method <- x$de_method[1]
    fgseaRes$reference_value <- x$reference_value[1]
    fgseaRes$coef_value <- project # Add grouping var to gsea df
    fgseaRes$formula_passed <- x$formula_passed[1]
    fgseaRes$formula <- x$formula[1]
    fgseaRes$cell_label_column <- x$cell_label_column[1]
    fgseaRes$cell_label <- x$cell_label[1]
    fgseaRes$sample_size <- sample_size
    fgseaRes$min_set_size <- min_set_size
    fgseaRes$max_set_size <- max_set_size
    fgseaRes$eps <- eps
    fgseaRes$score_type <- score_type
    fgseaRes$signed_ranking <- !arguments$options$unsigned_ranking

    ## Add categorical information for each id
    fgseaRes$gset_category_key <- gene_set_info[fgseaRes$annot_id,
                                                'gsetName']
    fgseaRes$gset_category <- gene_set_info[fgseaRes$annot_id,
                                            'gsetBioName']
    fgseaRes$gset_size <- Matrix::colSums(gene_set_genes[ , fgseaRes$annot_id])
    return(fgseaRes)
  })
gsea_rez <- do.call(rbind, fgsea_rl)
saveRDS(gsea_rez, file="fgsea.rds")

# Plot results
# vol_plot <- plot_volcano_plot(gsea_rez,
#                               'annot_coef',
#                               'pvalue',
#                               'gset_size',
#                               'gset_category',
#                               'coef_value',
#                               0.8,
#                               0.05,
#                               'annot_id')
# png(file = sprintf('%s-volcano_plot.png', output_file_base),
#     width = 12500,
#     height = 7500,
#     res = 600
# )
# print(vol_plot)
# dev.off()

# Plot results
# bubble_plot <- plot_bubble_plot(gsea_rez,
#                                 'annot_id',
#                                 'pvalue',
#                                 'gset_category',
#                                 'gset_size',
#                                 'coef_value')
# png(file = sprintf('%s-bubble_plot.png', output_file_base),
#     width = 12500,
#     height = 7500,
#     res = 600
# )
# print(bubble_plot)
# dev.off()

# if (fdr_n != 0) {
#   ## Only do FDR -- pvals are inflated
#   boxplot <- plot_boxplot(gsea_rez,
#                           'annot_coef',
#                           'annot_id',
#                           'fdr_pvalue_permute',
#                           'FDR',
#                           'coef_value',
#                           sig_threshold = 0.05)
#
#   png(file = sprintf('%s-boxplot.png', output_file_base),
#       units='px',
#       width=1600,
#       height=1600,
#       res=300
#   )
#   print(boxplot)
#   dev.off()
# }

# Save mappings
## First need to map `leadingEdge` to a string
gsea_rez$leadingEdge <- sapply(gsea_rez$leadingEdge, function(x) {
  return(paste(x, collapse = ","))
})

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

################################################################################
