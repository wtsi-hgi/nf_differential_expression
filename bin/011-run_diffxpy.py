#!/usr/bin/env python


__author__ = 'Monika Krzak and Leland Taylor'
__date__ = '2020-07-02'
__version__ = '0.0.1'

import argparse
import warnings
import random
from distutils.version import LooseVersion
import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import StandardScaler
import tensorflow as tf  # import before batchglm/diffxpy
import joblib  # for numpy matrix, joblib faster than pickle


# avoid tk.Tk issues
import matplotlib
matplotlib.use('Agg')


# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)
# 3. Set the `tensorflow` pseudo-random generator at a fixed value
if LooseVersion(tf.__version__) < '2.0.0':
    tf.random.set_random_seed(seed_value)
else:
    tf.random.set_seed(seed_value)


# Set up diffxpy backend
backend_diffxpy = 'numpy'

# Set up scaler for continuous_variables
scaler_continuous = StandardScaler(with_mean=True, with_std=True)


# Set GPU memory limits if running GPU
if LooseVersion(tf.__version__) <= '2.0.0':
    # with tf.Session() as sess:
    #     devices = sess.list_devices()
    # gpu = False
    # for dev in devices:
    #     if 'GPU' in dev.name:
    #         gpu = True
    #         break
    gpu = tf.test.is_gpu_available()
    if gpu:
        # If we can find gpu then use that as the diffxpy backend
        backend_diffxpy = 'tf1'

        # This enables tensorflow for diffxpy
        # More here:
        # https://diffxpy.readthedocs.io/en/latest/parallelization.html
        os.environ.setdefault('TF_NUM_THREADS', '1')
        os.environ.setdefault('TF_LOOP_PARALLEL_ITERATIONS', '1')

        # For TF v1
        # config = tf.ConfigProto()
        # config.gpu_options.allow_growth = True
        # session = tf.Session(config=config)
    else:
        warnings.warn('No GPUs detected.')
else:
    gpus = tf.config.list_physical_devices('GPU')
    print(gpus)
    if gpus:
        # If we can find gpu then use that as the diffxpy backend
        backend_diffxpy = 'tf2'

        # This enables tensorflow for diffxpy
        # More here:
        # https://diffxpy.readthedocs.io/en/latest/parallelization.html
        os.environ.setdefault('TF_NUM_THREADS', '1')
        os.environ.setdefault('TF_LOOP_PARALLEL_ITERATIONS', '1')

        # For TF v2
        try:
            # Method 1:
            # Currently, memory growth needs to be the same across GPUs
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)

            # Method 2:
            # Restrict TensorFlow to only allocate 1GB of memory on the first
            # GPU
            # tf.config.experimental.set_virtual_device_configuration(
            #     gpus[0],
            #     [tf.config.experimental.VirtualDeviceConfiguration(
            #         memory_limit=options.memory_limit*1024
            #     )])
            # logical_gpus = tf.config.list_logical_devices('GPU')
            # print(
            #     len(gpus),
            #     "Physical GPUs,",
            #     len(logical_gpus),
            #     "Logical GPUs"
            # )
        except RuntimeError as e:
            # Virtual devices must be set before GPUs have been initialized
            print(e)
    else:
        warnings.warn('No GPUs detected.')


import diffxpy.api as de  # import after tensorflow


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Differential expression analysis using Wald test.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-cond', '--condition_column',
        action='store',
        dest='condition_column',
        required=True,
        help='Condition column. Continuous variables are allowed, but this \
            has not been extensively tested.'
    )

    # parser.add_argument(
    #     '-cov', '--covariate_columns',
    #     action='store',
    #     dest='covariate_columns',
    #     default='',
    #     help='Comma seperated list of covariates (columns in adata.obs). \
    #         If this list includes normalization_factor then this covariate \
    #         is modeled as the size_factor. If normalization_factor is not \
    #         included then size factors are calculated as in \
    #         scater::librarySizeFactors. \
    #         Note: all continuous covariates are mean centered and \
    #         standardized. \
    #         (default: %(default)s)'
    # )

    parser.add_argument(
        '-cov_cont', '--covariate_columns_continuous',
        action='store',
        dest='covariate_columns_continuous',
        default='',
        help='Comma seperated list of covariates (columns in adata.obs). \
            If this list includes normalization_factor then this covariate \
            is modeled as the size_factor. If normalization_factor is not \
            included then size factors are calculated as in \
            scater::librarySizeFactors. \
            Note: all continuous covariates are mean centered and \
            standardized. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-cov_discrete', '--covariate_columns_discrete',
        action='store',
        dest='covariate_columns_discrete',
        default='',
        help='Comma seperated list of covariates (columns in adata.obs). \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-cl', '--cell_label_column',
        action='store',
        dest='cell_label_column',
        default='cluster',
        help='Anndata cell type label name in obs slot. (default: %(default)s)'
    )

    parser.add_argument(
        '-cla', '--cell_label_analyse',
        action='store',
        dest='cell_label_analyse',
        default='',
        help='Only analyse cells with these labels. Either "all" or comma \
            seperated list of cell labels. (default: %(default)s)'
    )

    parser.add_argument(
        '-m', '--method',
        action='store',
        dest='method',
        default='wald',
        help='Differential expression method to use. (default: %(default)s)'
    )

    # parser.add_argument(
    #     '-f', '--filter_low_genes',
    #     action='store',
    #     dest='filter_low_genes',
    #     default='wald',
    #     help='Filter for lowly expressed genes. (default: %(default)s)'
    # )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='diffxpy',
        help='Basename of output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the parameters
    condition_column = options.condition_column
    # covariate_columns_string = options.covariate_columns
    cell_label_column = options.cell_label_column
    cell_label_analyse = options.cell_label_analyse
    out_file_base = options.of

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Diffxpy automatically assigns covariates to categorical. List to stash
    # continuous variables.
    continuous_variables = []

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)
    adata.X = adata.layers['counts']

    # Sort out the condition we want to test
    adata.obs['condition'] = adata.obs[condition_column]
    # Check to see if the condition is a categorical variable
    if adata.obs[condition_column].dtype.name != 'category':
        # raise Exception('Condition is not a category.')
        continuous_variables.append('condition')
        warnings.warn('Treating condition as a continuous variable')
    # The below code converts the condition into an int and asserts it binary
    # adata.obs['condition'] = adata.obs[condition_column].cat.codes
    # if len(np.unique(adata.obs['condition']) != 2):
    #     raise Exception('There are not exactly 2 conditions.')

    # Figure out continuous covariates.
    # Also if covariate is also the condition, drop.
    covariate_columns = []
    # NOTE: below code is for if there is only one input covariate stream
    # if covariate_columns_string != '':
    #     for i in covariate_columns_string.split(','):
    #         if i == condition_column:
    #             msg = 'Condition ({}) is also covariate. {}.'.format(
    #                 i
    #             )
    #             raise Exception(msg)
    #             # warnings.warn(
    #             #     '{} {}.'.format(
    #             #         msg,
    #             #         'Dropping condition from the covariates list.'
    #             #     )
    #             # )
    #         else:
    #             covariate_columns.append(i)
    #             if adata.obs[i].dtype.name != 'category':
    #                 continuous_variables.append(i)
    # Add in continuous covariates
    if options.covariate_columns_continuous != '':
        for i in options.covariate_columns_continuous.split(','):
            # If pandas thinks a continuous variable is categorical, we have
            # a problem. Raise error.
            if adata.obs[i].dtype.name == 'category':
                raise Exception('Continuous variable appears to be category')
            covariate_columns.append(i)
            continuous_variables.append(i)
    # Add in discrete covariates
    if options.covariate_columns_discrete != '':
        for i in options.covariate_columns_discrete.split(','):
            # If a discrete columns int or float, cast to category.
            if adata.obs[i].dtype.name != 'category':
                adata.obs[i] = adata.obs[i].astype('category')
            covariate_columns.append(i)

    # Select the subset of cells we want for analysis
    if cell_label_analyse != 'all':
        cell_label_analyse = cell_label_analyse.split(',')
        adata = adata[
            adata.obs[cell_label_column].isin(cell_label_analyse)
        ].copy()  # NOTE: need copy here so it is not a view
    # clusters = np.sort(adata.obs[cell_label_column].unique())

    # Check to make sure that within each cluster, there are >1 values.
    if adata.obs[condition_column].dtype.name == 'category':
        n_cells_condition_cluster = adata.obs.groupby(
            [condition_column, cell_label_column]
        ).size().unstack()
        if n_cells_condition_cluster.values.min() == 0:
            raise Exception(
                'For one cell_label_column there are 0 conditions.'
            )
        if len(np.unique(adata.obs['condition'].cat.codes)) <= 1:
            raise Exception('There is only 1 condition.')

    # Check to make sure that for ecah covariate, there are >1 value.
    # NOTE: need a copy for loop because changing covariate_columns within
    # the loop.
    for i in covariate_columns.copy():
        unique_vals = adata.obs[i].unique().astype('str')
        # print('{}\t{}'.format(i, ''.join(unique_vals)))
        if len(unique_vals) <= 1:
            warnings.warn(
                'Removing covariate with only one value:\t{} {}'.format(
                    i,
                    unique_vals
                )
            )
            covariate_columns.remove(i)
            if adata.obs[i].dtype.name != 'category':
                continuous_variables.remove(i)

    # Get the size factors for analysis
    if 'normalization_factor' not in covariate_columns:
        # Calculate size factor covariate manually.
        #
        # When running diffxpy regress size_factors. For more details check:
        # https://nbviewer.jupyter.org/github/theislab/diffxpy_tutorials/blob/master
        # /diffxpy_tutorials/test/introduction_differential_testing.ipynb
        # "Inclusion of continuous effects"
        #
        # Option 1:
        # scater::library_size_factors
        adata.obs['diffxpy_size_factor'] = adata.obs[
            'total_counts'
        ] / adata.obs['total_counts'].mean()
        #
        # Option 2: ... same as above
        # https://github.com/theislab/diffxpy/issues/149
        # size_factors = np.asarray(adata.X.mean(axis=1))
        # size_factors = size_factors / size_factors.mean()
        # adata.obs['diffxpy_size_factor'] = size_factors
        # Option 3: Modified median of ratios method from Deseq2 from
        # # NOTE: this is not recommended since many of the final size factors
        # # will be 0.
        # # scater::medianSizeFactors
        # # See details here:
        # # https://rdrr.io/bioc/scater/man/medianSizeFactors.html
        # # https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
        # adata.X = adata.X.todense()
        # # Step 1: creates a pseudo-reference sample (geometric mean)
        # # For each gene a pseudo-reference sample is created that is equal to
        # # [original deseq] the geometric mean across all samples.
        # # pseudo_reference = sp.stats.gmean(adata.X, axis=0)
        # # [scater modification] the arithmetic mean.
        # pseudo_reference = adata.X.mean(axis=0)
        # # Step 2: calculates ratio of each sample to the reference
        # ratio_sample_reference = pd.DataFrame(adata.X / pseudo_reference)
        # # Step 3: calculate the normalization factor for each sample
        # # (size factor)
        # # if np.any(np.all(np.isnan(ratio_sample_reference), axis=1)):
        # #     raise Exception('Error in size factor calculation')
        # size_factor = ratio_sample_reference.median(axis=1)
        # # Step 4: calculate the normalized count values using the
        # # normalization factor
        # # normalized_counts = adata.X / size_factor
    else:
        adata.obs['diffxpy_size_factor'] = adata.obs['normalization_factor']
        continuous_variables.remove('normalization_factor')
        covariate_columns.remove('normalization_factor')

    # TODO: Check for nan in covariates and conditions?

    # Center and standardize continuous variables
    if len(continuous_variables) > 0:
        adata.obs[continuous_variables] = scaler_continuous.fit_transform(
            adata.obs[continuous_variables]
        )
        print(adata.obs[continuous_variables].mean(axis=0))
        print(adata.obs[continuous_variables].std(axis=0))
        print('Continuous varibles:\t{}'.format(
            ','.join(continuous_variables))
        )

    # TODO: filter lowly expressed genes for this cluster?
    # Lowly expressed genes tend to result in false postives because outliers
    # pull up the regression line. These, however, can always be filtered out
    # afterwards.

    # Describe the input covariates
    if len(covariate_columns) > 0:
        print('Covariate description:\t{}'.format(','.join(covariate_columns)))
        print(adata.obs[covariate_columns].describe(include='all'))

    # NOTE:
    # One may get a constrained design matrix is not full rank when for
    # the categorcial variable combinations, one facet is all ones.
    discrete_vars = [
        i for i in covariate_columns if i not in continuous_variables
    ]
    if len(discrete_vars) > 1:
        n_cells_condition = adata.obs.groupby(discrete_vars).size().unstack()
        print(n_cells_condition)

    # Run diffxpy
    #
    # Tutorial for conditions:
    # https://nbviewer.jupyter.org/github/theislab/diffxpy_tutorials/blob/master/diffxpy_tutorials/test/introduction_differential_testing.ipynb
    # Tutorial for continuous:
    # https://nbviewer.jupyter.org/github/theislab/diffxpy_tutorials/blob/master/diffxpy_tutorials/test/modelling_continuous_covariates.ipynb
    formula = '~ 1 + condition'
    for cov in covariate_columns:
        formula = '{} + {}'.format(formula, cov)
    print('Formula:\t{}'.format(formula))
    # Get the coefficient names.
    # Example if coefficient = M or F then...
    # coef_to_test=['condition[T.M]']
    # If T1, T2, T3
    # coef_to_test=['condition[T.1]', 'condition[T.2]', 'condition[T.3]']
    coef_names = de.utils.preview_coef_names(
        sample_description=adata.obs,
        formula=formula,
        as_numeric=continuous_variables
    )
    coef_names = [i for i in coef_names if i.startswith('condition')]
    # param as_numeric:
    # Which columns of sample_description to treat as numeric and
    # not as categorical. This yields columns in the design matrix
    # which do not correspond to one-hot encoded discrete factors.
    # This makes sense for number of genes, time, pseudotime or space
    # for example.
    print('Suggested backend:\t{}'.format(backend_diffxpy))
    warnings.warn('Backend hardset to numpy due to bug in diffxpy.')
    print('Using backend:\t{}'.format('numpy'))
    if options.method == 'wald':
        test_results = de.test.wald(
            data=adata,
            formula_loc=formula,
            factor_loc_totest=coef_names,
            as_numeric=continuous_variables,
            size_factors='diffxpy_size_factor',
            backend='numpy',  # numpy, tf1, or tf2
            noise_model='nb',
            # batch_size=100,  # the memory load of the fitting procedure
            quick_scale=False,
        )
    # elif options.method == 'lrt':
    #     test_results = de.test.lrt(
    #         data=adata,
    #         full_formula_loc=formula,
    #         reduced_formula_loc=formula.replace('+ condition', ''),
    #         # factor_loc_totest=coef_names,
    #         as_numeric=continuous_variables,
    #         size_factors='diffxpy_size_factor',
    #         backend=backend_diffxpy,  # numpy, tf1, or tf2
    #         noise_model='nb',
    #         # batch_size=100,  # the memory load of the fitting procedure
    #         quick_scale=False
    #     )
    else:
        raise Exception('Invalid method.')

    # Make a table of the results
    df_results = test_results.summary()
    df_results = df_results.rename(columns={'coef_mle': 'coef'})
    df_results['de_method'] = 'diffxpy-{}'.format(options.method)
    df_results['condition'] = condition_column
    coef_names.sort()
    df_results['coef_condition'] = ','.join(coef_names)
    # covariates_passed = covariates that were passed to this script call.
    covariates_passed = options.covariate_columns_discrete.split(',')
    if options.covariate_columns_continuous != '':
        for i in options.covariate_columns_continuous.split(','):
            covariates_passed.append(i)
    covariates_passed.sort()
    df_results['covariates_passed'] = ','.join(covariates_passed)
    # covariates = covariates that were actually included in the model since
    # covariates with only one value are dropped.
    covariate_columns.sort()
    df_results['covariates'] = ','.join(covariate_columns)
    df_results['cell_label_column'] = cell_label_column
    cell_label_analyse.sort()
    df_results['cell_label_analysed'] = ','.join(cell_label_analyse)
    df_results = df_results.sort_values(
        by=['pval', 'log2fc', 'mean'],
        ascending=[True, False, False]
    )
    df_results.to_csv(
        '{}-de_results.tsv.gz'.format(out_file_base),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )

    # Now save the results object
    out_f = '{}-de_results_obj.joblib.gz'.format(out_file_base)
    joblib.dump(
        test_results,
        out_f,
        compress=('gzip', 3)
    )

    # Make plots
    _ = test_results.plot_ma(
        corrected_pval=False,
        show=False,
        log2_fc_threshold=10,
        min_mean=0.0001,
        alpha=0.05,
        size=20,
        save='{}'.format(out_file_base),
        suffix='-plot_ma.png',
        return_axs=False
    )
    test_results.plot_volcano(
        corrected_pval=False,
        show=False,
        log10_p_threshold=-30,
        log2_fc_threshold=10,
        alpha=0.05,
        min_fc=1,
        size=20,
        save='{}'.format(out_file_base),
        suffix='-plot_volcano.png',
        return_axs=False
    )

    # NOTE: the below code partitions the data within diffxpy ... we don't
    # need to do this since we have already partitioned clusters via the
    # cell_label_analyse command
    # part = de.test.partition(data=adata, parts="cluster")
    #
    # if options.cov == "":
    #     formula = "~ 1 + condition"
    #     test_part = part.wald(formula_loc=formula,
    #                           factor_loc_totest="condition")
    # else:
    #     formula = "~ 1 + condition + " + options.cov
    #     test_part = part.wald(formula_loc=formula,
    #                           factor_loc_totest="condition",
    #                           as_numeric=[options.cov])
    # for i in clusters:
    #     df = test_part.tests[test_part.partitions.index(i)].summary()
    #     df.to_csv('{}-cluster_de_wald.tsv'.format(i), sep='\t', index=False)
    #
    # part = de.test.partition(data=adata, parts="cluster")
    # test_part = part.wald(formula_loc="~ 1 + condition",
    #                       factor_loc_totest="condition")


if __name__ == '__main__':
    main()
