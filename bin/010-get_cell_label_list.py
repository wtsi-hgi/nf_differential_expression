#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-07-14'
__version__ = '0.0.1'

import argparse
import numpy as np
import scanpy as sc


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Get a list of all celltypes in an anndataframe.
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
        '--cell_label',
        action='store',
        dest='cell_label',
        default='cluster',
        help='Anndata cell type label name in obs slot. (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    np_array = np.sort(adata.obs[options.cell_label].unique().astype(str))
    np.savetxt('cell_labels.csv', np_array, fmt='%s')


if __name__ == '__main__':
    main()
