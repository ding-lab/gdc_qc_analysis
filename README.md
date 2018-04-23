```
python make_db.py --db-url 'sqlite:///all_variants.sqlite' \
    --mc3-maf ../GDC_QC_data/MC3/mc3.v0.2.8.PUBLIC.GRCh38_converted.maf.gz \
    --gdc-root ../GDC_QC_data/GDC_data_release/Release_10.0
/usr/local/Cellar/sqlite/3.23.1/bin/sqlite3 all_variants.sqlite < group_gdc_callers.sql
/usr/local/Cellar/sqlite/3.23.1/bin/sqlite3 all_variants.sqlite < mc3_select_samples.sql
```
