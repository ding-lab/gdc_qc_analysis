import argparse
import logging
from pathlib import Path
from sqlalchemy import (
    create_engine, event,
    MetaData, Table, Column, Integer, Text, Index,
    UniqueConstraint
)
from sqlalchemy.engine import Engine
from maf_utils import GDCMAF, MC3MAF

logger = logging.getLogger(__name__)
BATCH_SIZE = 1000


def define_db_schema(metadata, mc3_maf, gdc_maf):
    mc3_cols_integer = [
        'start_position',
        'end_position',
        'strand_vep',
        'raw_file_line_number',
    ]
    mc3_cols = []
    for col in mc3_maf.columns:
        if col in mc3_cols_integer:
            c = Column(col, Integer())
        else:
            c = Column(col, Text())
        mc3_cols.append(c)
    Table(
        'mc3', metadata,
        *mc3_cols,
        Index('mc3_ix_tumor_barcode', 'tumor_sample_barcode'),
        # Unique constraint
        UniqueConstraint('raw_file_line_number'),
    )

    gdc_cols_integer = [
        'start_position',
        'end_position',
        'raw_file_line_number',
    ]
    gdc_cols = []
    for col in gdc_maf.columns:
        if col in gdc_cols_integer:
            c = Column(col, Integer())
        else:
            c = Column(col, Text())
        gdc_cols.append(c)
    Table(
        'gdc', metadata,
        *gdc_cols,
        Index('gdc_ix_tumor_barcode', 'tumor_sample_barcode'),
        # Unique constraint
        UniqueConstraint('cancer_type', 'caller', 'raw_file_line_number'),
    )


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA cache_size=-4192000")
    cursor.execute("PRAGMA temp_store=MEMORY")
    cursor.execute("PRAGMA journal_mode=MEMORY")
    cursor.close()


def load_maf(conn, metadata, maf, db_table):
    ins = db_table.insert()
    ins_batch = []
    for i, record in enumerate(maf, 1):
        if i % 100000 == 0:
            logger.info(f'... inserted {i:,d} records')
        if len(ins_batch) >= BATCH_SIZE:
            with conn.begin():
                conn.execute(ins, ins_batch)
            ins_batch = []
        ins_batch.append(record._asdict())

    # Load the last batch less than the batch size
    if ins_batch:
        conn.execute(ins, ins_batch)


def main(db_url, mc3_maf_pth, gdc_root):
    # Read all MAFs
    mc3_maf = MC3MAF(Path(mc3_maf_pth))
    gdc_maf_pths = list(Path(gdc_root).glob('*/TCGA.*.somatic.maf.gz'))
    gdc_mafs = [GDCMAF(pth) for pth in gdc_maf_pths]

    # Create database schema
    metadata = MetaData()
    db_engine = create_engine(db_url)

    define_db_schema(metadata, mc3_maf, gdc_mafs[0])
    metadata.create_all(db_engine, checkfirst=True)

    # Load in data
    logger.info(f'Load variants to {db_url}')
    conn = db_engine.connect()

    logger.info(f'Loading MC3 variants')
    load_maf(conn, metadata, mc3_maf, metadata.tables['mc3'])

    logger.info(f'Loading GDC variants')
    for maf in gdc_mafs:
        logger.info(f'Loading GDC {maf.cancer_type} {maf.caller}')
        load_maf(conn, metadata, maf, metadata.tables['gdc'])

    logger.info(f'All variants are loaded to {db_url}')


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
        '--db-url', required=True,
        default='sqlite:///variants.sqlite',
        help="Database connection URL"
    )
    parser.add_argument(
        '--mc3-maf', required=True,
        help='Path to the hg38 version of MC3 MAF'
    )
    parser.add_argument(
        '--gdc-root', required=True,
        help='Path to the GDC data release root folder'
    )
    return parser


if __name__ == '__main__':
    parser = setup_cli()
    args = parser.parse_args()

    main(args.db_url, args.mc3_maf, args.gdc_root)
