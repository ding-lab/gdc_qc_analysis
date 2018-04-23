Conda environment setup:

    conda create -n gdc_qc python=3.6 notebook matplotlib sqlalchmey pandas

The script to generate the database:

```
python make_db.py --db-url 'sqlite:///processed_data/all_variants.sqlite' \
    --mc3-maf ../GDC_QC_data/MC3/mc3.v0.2.8.PUBLIC.GRCh38_converted.maf.gz \
    --gdc-root ../GDC_QC_data/GDC_data_release/Release_10.0
sqlite3 all_variants.sqlite < group_gdc_callers.sql
sqlite3 all_variants.sqlite < mc3_select_samples.sql
```
