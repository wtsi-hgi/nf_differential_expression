#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-08-19'
__version__ = '0.0.1'

import argparse
import scanpy as sc
import plotnine as plt9
import pandas as pd
import os
from scipy import stats
import numpy as np
import random


# Set seed for reproducibility
seed_value = 0
os.environ['PYTHONHASHSEED'] = str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Variance decomposition.
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
        '-cov', '--covariates',
        action='store',
        dest='cov',
        required=True,
        help='Covariates to be decomposed.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='variance_decomposition',
        help='Basename of output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    out_file_base = options.of

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    # sc.pp.subsample(adata, n_obs=5000, random_state=0, copy=False)

    df = adata.obs

    covariates = options.cov.split(',')
    covariates = np.asarray(covariates)

    # covariates = ['sex', 'age', 'smoking_status', 'disease_status',
    #               'experimentalist', 'sanger_sample_id', 'cluster',
    #               'total_counts', 'total_genes_detected', 'id_run',
    #               'lane', 'bead_version', 'bead_lot', 'chip_version',
    #               'chip_lot']

    for j in covariates:
        if df[j].dtype.name == 'category':
            df[j] = pd.get_dummies(data=df[j], drop_first=True)

    # Run variance decomposition
    M = adata.layers['log1p_cp10k'].todense().copy()

    result = []
    for j in covariates:
        print(j)

        y = np.array([])

        for i in range(0, M.shape[1]):

            slope, intercept, r_value, p_value, std_err = stats.linregress(
                                    np.array(M[:, i].ravel()), np.array(df[j]))
            y = np.append(y, r_value**2)

        # cor_value = stats.spearmanr(M[:,i], df[j])
        # y = np.append(y, cor_value)

        result.append(y)

    df1 = pd.DataFrame(np.transpose(result), columns=covariates)
    df1.to_csv(
        '{}-matrix.tsv.gz'.format(out_file_base),
        sep='\t',
        compression='gzip',
        index=False)

    df1 = pd.melt(df1)
    df1['value'] = df1['value']*100

    # Density plot  - all in one
    gplt = plt9.ggplot(df1, plt9.aes(x='value', color='variable'))
    gplt = gplt + plt9.geom_line(stat='density', alpha=0.7, size=1)
    gplt = gplt + plt9.scale_x_log10(
        breaks=np.array([0.001, 0.01, 0.1, 1, 10, 100]),
        labels=np.array([0.001, 0.01, 0.1, 1, 10, 100]),
        limits=np.array([0.001, 100]))
    gplt = gplt + plt9.ylab("density") + plt9.xlab("% of variance explained")
    plt9.ggsave(plot=gplt, filename="variance_decomposed_density_plot",
                width=12, height=8)

    # Density plot  - facet by covariates
    gplt1 = plt9.ggplot(df1, plt9.aes(x='value', color='variable'))
    gplt1 = gplt1 + plt9.geom_line(stat='density', alpha=0.7, size=1)
    gplt1 = gplt1 + plt9.facet_wrap('variable', scales='free_x')
    gplt1 = gplt1 + plt9.theme_bw()
    gplt1 = gplt1 + plt9.scale_x_log10(
        breaks=np.array([0.001, 0.01, 0.1, 1, 10, 100]),
        labels=np.array([0.001, 0.01, 0.1, 1, 10, 100]),
        limits=np.array([0.001, 100]))
    gplt1 = gplt1 + plt9.ylab("density") + plt9.xlab("% of variance explained")
    plt9.ggsave(plot=gplt1, filename="variance_decomposed_density_plot_facet",
                width=12, height=8)

    # Boxplot
    gplt2 = plt9.ggplot(df1, plt9.aes(x='variable', y='value', color='variable')
                        )
    gplt2 = gplt2 + plt9.geom_violin()
    gplt2 = gplt2 + plt9.theme_bw()
    gplt2 = gplt2 + plt9.scale_x_log10(
        breaks=np.array([0.001, 0.01, 0.1, 1, 10, 100]),
        labels=np.array([0.001, 0.01, 0.1, 1, 10, 100]),
        limits=np.array([0.001, 100]))
    gplt2 = gplt2 + plt9.ylab("density") + plt9.xlab("% of variance explained")
    gplt2 = gplt2 + plt9.geom_boxplot(width=0.1)
    plt9.ggsave(plot=gplt2, filename="variance_decomposed_boxplot",
                width=12, height=8)


if __name__ == '__main__':
    main()
