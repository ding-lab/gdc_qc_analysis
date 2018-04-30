configfile: 'config.yaml'

from pathlib import Path

CANCER_TYPES = ['BRCA', 'COAD', 'LAML', 'OV']
GDC_CALLERS = ['mutect', 'somaticsniper', 'muse', 'varscan']


def find_all_gdc_mafs(gdc_root: Path):
    """Check if all the GDC MAFs exist."""
    gdc_mafs = []
    for cancer in CANCER_TYPES:
        for caller in GDC_CALLERS:
            try:
                maf = next(gdc_root.glob(f'*/TCGA.{cancer}.{caller}.*.somatic.maf.gz'))
            except StopIteration:
                raise ValueError(f'Cannot find GDC MAF for {cancer} by {caller}')
            gdc_mafs.append(maf)
    return gdc_mafs

GDC_MAFS = find_all_gdc_mafs(Path(config['GDC_DATA_ROOT']))


def find_mc3_maf(wildcards):
    """Return different MC3 MAFS

    Now it can be different access type (public or control).
    """
    return config['MC3_MAF_PTHS'][wildcards.access_type]


rule gen_g_coord_bed:
    """Generate coordinates only BED file."""
    input: maf=find_mc3_maf
    output: 'processed_data/mc3.{access_type}.GRCh37.g_coords.gz'
    shell:
        'python scripts/gen_g_coord_bed.py {input.maf} {output}'


rule g_coords_b37_to_b38:
    """Convert coordinates from GRCh37 to GRCh38."""
    input: 'processed_data/{name}.GRCh37.g_coords.gz'
    output: 'processed_data/{name}.GRCh38.g_coords.gz'
    params:
        crossmap_bin=config['CROSS_MAP_BIN'],
        chain_file=config['CHAIN_PTH']
    shell:
        '''
        {params.crossmap_bin} bed {params.chain_file} {input} \
            | gzip -c > {output}
        '''


rule swap_maf_g_coords:
    """Replace the MAF file with GRCh38 coordinates."""
    input:
        maf=find_mc3_maf,
        g_coords='processed_data/mc3.{access_type}.GRCh38.g_coords.gz'
    output: 'processed_data/mc3.{access_type}.converted.GRCh38.maf.gz'
    shell:
        'python scripts/swap_g_coord_bed.py {input.maf} {input.g_coords} {output}'


rule make_db:
    """Generate the SQLite database."""
    input:
        mc3_maf='processed_data/mc3.public.converted.GRCh38.maf.gz',
        gdc_mafs=GDC_MAFS
    params:
        gdc_root=config['GDC_DATA_ROOT']
    output: 'processed_data/variants.sqlite'
    run:
        shell("python scripts/make_db.py --db-url 'sqlite:///{output}' --mc3-maf {input.mc3_maf} --gdc-root {params.gdc_root}")
        shell('sqlite {output} < scripts/group_gdc_callers.sql')

        import sqlite3
        conn = sqlite3.connect(output[0])
        sql_cmd = 'SELECT DISTINCT tumor_sample_barcode FROM gdc_grouped_callers'
        tumor_barcodes = sorted(t[0] for t in conn.execute(sql_cmd).fetchall())

        TPL = '''\
        DROP TABLE IF EXISTS mc3_selected;
        CREATE TABLE IF NOT EXISTS mc3_selected AS
        SELECT * FROM mc3
        WHERE tumor_sample_barcode IN ({sample_list});
        DROP TABLE mc3;
        VACUUM;
        '''
        with open('scripts/mc3_select_samples.sql', 'w') as f:
            sample_list = ', '.join(f"'{v}'" for v in tumor_barcodes)
            print(TPL.format(sample_list=sample_list), file=f)

        shell('sqlite {output} < scripts/mc3_select_samples.sql')

rule all:
    input:
        'processed_data/mc3.public.converted.GRCh38.maf.gz',
        # 'processed_data/mc3.controlled.converted.GRCh38.maf.gz'
