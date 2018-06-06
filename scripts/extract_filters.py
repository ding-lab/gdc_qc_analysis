import argparse
from collections import Counter
from functools import partial
from itertools import chain
import logging

import pandas as pd
import numpy as np


logger = logging.getLogger(__name__)


def filter_val_indictor(filter_vals, all_possible_vals):
    """Convert the filter value into the order of all_possible_vals

    >>> filter_list_to_indicator(['B', 'C'], ['A', 'B', 'C'])
    [0, 1, 1]
    """
    return [val in filter_vals for val in all_possible_vals]


def gen_filter_indicator_df(df, filter_col, filter_name, sep=';'):
    """
    Convert a dense filter column into indicator columns ordered by the
    occurence frequency.


    For example, a given filter column in the data frame::

        filter_col
        ----------
        a;b;c
        c

    the output will become::

        filter_name__c  filter_name__b  filter_name__a
        --------------  --------------  --------------
        1               1               1
        1               0               0

    """
    # Replace NA will empty string and split the compressed filter
    splitted_filter_vals = [x.split(sep) for x in df[filter_col].fillna(value='')]

    # Count the filter occurence and use the freq to sort filter
    filter_counter = Counter(chain.from_iterable(splitted_filter_vals))

    # Collect all the filters and set their column names
    ordered_filter_vals = [kv_pair[0] for kv_pair in filter_counter.most_common()]
    ordered_filter_reprs = [f'{filter_name}__{val}' for val in ordered_filter_vals]

    # Convert the filter to indicator columns
    filter_indicator_df = pd.DataFrame(
        list(map(
            partial(filter_val_indictor, all_possible_vals=ordered_filter_vals),
            splitted_filter_vals
        )),
        columns=ordered_filter_reprs
    ).astype(int)
    return filter_indicator_df


def setup_cli():
    # Setup console logging
    console = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(console)
    log_fmt = '[%(asctime)s][%(levelname)-7s] %(message)s'
    log_formatter = logging.Formatter(log_fmt, '%Y-%m-%d %H:%M:%S')
    console.setFormatter(log_formatter)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('tsv', help="Path to the variant TSV")
    parser.add_argument('out', help="Output path")
    return parser


def main(tsv_pth, out_pth):
    df = pd.read_table(tsv_pth)
    gdc_filter_df = gen_filter_indicator_df(df, 'gdc_filter', 'gdc')
    gdc_gdc_filter_df = gen_filter_indicator_df(df, 'gdc_gdc_filter', 'gdc_gdc')
    mc3_filter_df = gen_filter_indicator_df(df, 'mc3_filter', 'mc3', sep=',')
    all_filter_df = pd.concat([gdc_filter_df, gdc_gdc_filter_df, mc3_filter_df], axis='columns')

    all_filter_df.to_csv(out_pth, sep='\t', index=False)


if __name__ == '__main__':
    parser = setup_cli()
    args = parser.parse_args()

    main(args.tsv, args.out)
