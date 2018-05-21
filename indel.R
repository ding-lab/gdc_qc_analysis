library(tidyverse)
library(DBI)
library(RSQLite)

conn = dbConnect(RSQLite::SQLite(), dbname='./processed_data/all_variants.somatic_only.sqlite')
dbExecute(conn, 'PRAGMA cache_size=-8000000')
dbExecute(conn, 'PRAGMA temp_store=MEMORY')
overlap_tbl <- as.tibble(dbReadTable(conn, 'full_overlap'))


indel_tbl <- overlap_tbl %>% filter(
    mc3_variant_type %in% c('DEL', 'INS') | gdc_variant_type %in% c('DEL', 'INS')
)

