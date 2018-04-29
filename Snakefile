from pathlib import Path


MC3_MAF_PTHS = {
    'public': '../GDC_QC_data/MC3/mc3.v0.2.8.PUBLIC.maf.gz',
    'controlled': '../GDC_QC_data/MC3/mc3.v0.2.8.CONTROLLED.CT.maf.gz'
}
CROSS_MAP_BIN = '~/miniconda3/envs/crossmap/bin/CrossMap.py'
CHAIN_PTH = '../GDC_QC_data/liftOver_chains/GRCh37_to_GRCh38.chain.gz'


rule gen_g_coord_bed:
    """Generate coordinate only BED file."""
    input: maf=lambda wildcards: MC3_MAF_PTHS[wildcards.access_type]
    output: 'processed_data/mc3.{access_type}.GRCh37.g_coords.gz'
    shell:
        'python scripts/gen_g_coord_bed.py {input.maf} {output}'


rule g_coords_b37_to_b38:
    input: 'processed_data/{name}.GRCh37.g_coords.gz'
    output: 'processed_data/{name}.GRCh38.g_coords.gz'
    params:
        crossmap_bin=CROSS_MAP_BIN,
        chain_file=CHAIN_PTH
    shell:
        '''
        {params.crossmap_bin} bed {params.chain_file} {input} \
            | gzip -c > {output}
        '''

rule swap_maf_g_coords:
    input:
        maf=lambda wildcards: MC3_MAF_PTHS[wildcards.access_type],
        g_coords='processed_data/mc3.{access_type}.GRCh38.g_coords.gz'
    output: 'processed_data/mc3.{access_type}.converted.GRCh38.maf.gz'
    shell:
        'python scripts/swap_g_coord_bed.py {input.maf} {input.g_coords} {output}'

rule all:
    input:
        'processed_data/mc3.public.converted.GRCh38.maf.gz',
        'processed_data/mc3.controlled.converted.GRCh38.maf.gz'
