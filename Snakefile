configfile: 'config.yaml'

from pathlib import Path

CANCER_TYPES = ['BRCA', 'COAD', 'LAML', 'OV']
GDC_CALLERS = ['mutect', 'somaticsniper', 'muse', 'varscan']


def find_all_gdc_mafs(gdc_root: Path, file_type='somatic'):
    """Check if all the GDC MAFs exist."""
    gdc_mafs = []
    for cancer in CANCER_TYPES:
        for caller in GDC_CALLERS:
            try:
                maf = next(gdc_root.glob(f'*/TCGA.{cancer}.{caller}.*.{file_type}.maf.gz'))
            except StopIteration:
                raise ValueError(f'Cannot find GDC MAF for {cancer} by {caller}')
            gdc_mafs.append(str(maf))
    return gdc_mafs

GDC_MAFS = find_all_gdc_mafs(Path(config['GDC_DATA_ROOT']))
GDC_PROTECTED_MAFS = find_all_gdc_mafs(
    Path(config['GDC_DATA_ROOT']), 'protected')


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
        'python scripts/swap_g_coord_bed.py {input.maf} {input.g_coords} | gzip -c > {output}'


rule make_db:
    """Generate the SQLite database."""
    input:
        mc3_maf='processed_data/mc3.public.converted.GRCh38.maf.gz',
        gdc_mafs=GDC_MAFS
    params:
        gdc_root=config['GDC_DATA_ROOT']
    output: 'processed_data/all_variants.sqlite'
    run:
        shell("python scripts/make_db.py --db-url 'sqlite:///{output}' --mc3-maf {input.mc3_maf} --gdc-root {params.gdc_root}")
        shell('sqlite3 -echo {output} < scripts/group_gdc_callers.sql')
        shell('sqlite3 -echo {output} < scripts/subset_samples.sql')
        shell("python scripts/create_overlap_table.py --db-pth {output}")
        # shell('sqlite3 -echo {output} < scripts/clean_up.sql')


rule add_protected_mafs_to_db:
    """Add protected MAFs into SQLite database."""
    input:
        mc3_maf='processed_data/mc3.controlled.converted.GRCh38.maf.gz',
        gdc_mafs=GDC_MAFS,
        db=rules.make_db.output[0]
    output:
        'processed_data/db_state/has_added_protected_mafs'
    params:
        gdc_root=config['GDC_DATA_ROOT']
    run:
        shell("python scripts/add_protected_maf.py --db-url 'sqlite:///{input.db}' --mc3-maf {input.mc3_maf} --gdc-root {params.gdc_root}")
        shell("sqlite3 -echo {input.db} < scripts/group_protected_gdc_callers_loose.sql")
        shell("sqlite3 -echo {input.db} < scripts/create_recoverable_unique_tables.sql")
        shell("touch {output}")


rule extract_recoverable_tables:
    input:
        db='processed_data/all_variants.sqlite',
        db_state='processed_data/db_state/has_added_protected_mafs'
    output:
        expand('processed_data/{grp}_recoverable_unique_variants.tsv.gz', grp=['gdc', 'mc3']),
        expand('processed_data/{grp}_not_recoverable_unique_variants.tsv.gz', grp=['gdc', 'mc3']),
    shell:
        'sqlite3 {input.db} < scripts/extract_recoverable_tables.sql'



rule extract_filter:
    input: 'processed_data/{grp}_recoverable_unique_variants.tsv.gz'
    output: 'processed_data/{grp}_recoverable_unique_variants.filter_cols.tsv.gz'
    shell:
        'python scripts/extract_filters.py {input} {output}'


rule all:
    input:
        'processed_data/mc3.public.converted.GRCh38.maf.gz',
        # 'processed_data/mc3.controlled.converted.GRCh38.maf.gz',
        'processed_data/db_state/has_added_protected_mafs',
        expand('processed_data/{grp}_recoverable_unique_variants.filter_cols.tsv.gz', grp=['gdc', 'mc3']),
