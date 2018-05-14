import argparse
import logging
import sqlite3


logger = logging.getLogger(__name__)

FIELDS_TO_REPORT = '''\
    chromosome, start_position, end_position, tumor_sample_barcode,
    reference_allele,
    g.tumor_seq_allele1 AS gdc_tumor_seq_allele1,
    m.tumor_seq_allele1 AS mc3_tumor_seq_allele1,
    tumor_seq_allele2,
    g.hgvsp AS gdc_hgvsp,
    m.hgvsp AS mc3_hgvsp,
    g.hgvsc AS gdc_hgvsc,
    m.hgvsc AS gdc_hgvsc,
    (CASE WHEN g.exon='.' THEN NULL ELSE g.exon END) AS gdc_exon,
    (CASE WHEN m.exon='.' THEN NULL ELSE m.exon END) AS mc3_exon,
    g.variant_classification AS gdc_variant_classification,
    m.variant_classification AS mc3_variant_classification,
    g.variant_type AS gdc_variant_type,
    m.variant_type AS mc3_variant_type,
    g.hugo_symbol AS gdc_hugo_symbol,
    m.hugo_symbol AS mc3_hugo_symbol,
    m.transcript_id AS mc3_transcript_id,
    g.transcript_id AS gdc_transcript_id,
    g.callers AS gdc_callers,
    m.centers AS mc3_callers,
    m.ncallers AS mc3_ncallers,
    m.filter AS mc3_filter,
    g.filter AS gdc_filter,
    g.gdc_filter AS gdc_gdc_filter,
    g.mc3_overlap AS gdc_mc3_overlap,
    g.gdc_validation_status AS gdc_validation_status,
    g.t_depth_per_caller AS gdc_t_depth_per_caller,
    g.t_ref_count_per_caller AS gdc_t_ref_count_per_caller,
    g.t_alt_count_per_caller AS gdc_t_alt_count_per_caller,
    g.n_depth_per_caller AS gdc_n_depth_per_caller,
    g.n_ref_count_per_caller AS gdc_n_ref_count_per_caller,
    g.n_alt_count_per_caller AS gdc_n_alt_count_per_caller,
    m.t_depth AS mc3_t_depth, m.t_ref_count AS mc3_t_ref_count, m.t_alt_count AS mc3_t_alt_count,
    m.n_depth AS mc3_n_depth, m.n_ref_count AS mc3_n_ref_count, m.n_alt_count AS mc3_n_alt_count,
    g.context AS gdc_context,
    m.context AS mc3_context,
    (CASE WHEN g.existing_variation='' THEN NULL ELSE g.existing_variation END) AS gdc_existing_variation,
    (CASE WHEN m.existing_variation='.' THEN NULL ELSE m.existing_variation END) AS mc3_existing_variation,
    m.rowid AS mc3_rowid,
    g.rowid AS gdc_rowid,
    (CASE WHEN m.rowid IS NOT NULL AND g.rowid IS NOT NULL THEN 1 ELSE 0 END) AS shared_by_gdc_mc3,
    (CASE WHEN m.rowid IS NULL AND g.rowid IS NOT NULL THEN 1 ELSE 0 END) AS only_in_gdc,
    (CASE WHEN m.rowid IS NOT NULL AND g.rowid IS NULL THEN 1 ELSE 0 END) AS only_in_mc3
'''


def main(db_pth):
    conn = sqlite3.connect(db_pth)
    conn.executescript('''\
    PRAGMA cache_size=-4192000;
    PRAGMA temp_store=MEMORY;
    PRAGMA journal_mode=OFF;
    ''')

    conn.executescript(f'''\
    DROP TABLE IF EXISTS full_overlap;

    CREATE TABLE IF NOT EXISTS full_overlap AS
    WITH sample_cancer_type AS (
        SELECT DISTINCT tumor_sample_barcode, cancer_type
        FROM gdc_shared_samples
    ), full_overlap AS (
        SELECT {FIELDS_TO_REPORT}
        FROM gdc_shared_samples g
        LEFT JOIN mc3_shared_samples m
            USING (tumor_sample_barcode, chromosome, start_position, end_position, reference_allele, tumor_seq_allele2)
        UNION ALL
        SELECT {FIELDS_TO_REPORT}
        FROM mc3_shared_samples m
        LEFT JOIN gdc_shared_samples g
            USING (tumor_sample_barcode, chromosome, start_position, end_position, reference_allele, tumor_seq_allele2)
        WHERE g.tumor_seq_allele1 IS NULL
    )
    SELECT full_overlap.*, sample_cancer_type.cancer_type
    FROM full_overlap
    LEFT JOIN sample_cancer_type
        USING (tumor_sample_barcode)
    ;

    CREATE INDEX ix_full_overlap_tumor_sample_barcode ON full_overlap (tumor_sample_barcode);
    '''
    # CREATE INDEX ix_full_overlap_genom_range ON full_overlap (
    #     chromosome, start_position, end_position DESC
    # );
    )


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
    parser.add_argument(
        '--db-pth', required=True,
        default='./variants.sqlite',
        help="Path to the SQLite database"
    )
    return parser


if __name__ == '__main__':
    parser = setup_cli()
    args = parser.parse_args()

    main(args.db_pth)
