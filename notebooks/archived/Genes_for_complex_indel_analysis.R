library(tidyverse)
library(DBI)
library(RSQLite)

conn = dbConnect(RSQLite::SQLite(), dbname='./processed_data/all_variants.sqlite')
dbExecute(conn, 'PRAGMA cache_size=-4192000')
dbExecute(conn, 'PRAGMA temp_store=MEMORY')

tp53_tbl <- as.tibble(dbGetQuery(
    conn, 
    'SELECT * FROM full_overlap WHERE mc3_hugo_symbol = \'TP53\' OR gdc_hugo_symbol = \'TP53\''
))

write_tsv(tp53_tbl, '~/Box Sync/Ding_Lab/Projects_Current/GDC-QC/Results/TP53_overlap.tsv')
