PRAGMA cache_size=-8000000;  -- Use 8GB RAM as cache
PRAGMA temp_store=MEMORY;
.header on
.mode csv
.separator \t

-- GDC recoverable unique
.once '| gzip -c > processed_data/gdc_recoverable_unique_variants.tsv.gz'
SELECT g.*,
    m.hugo_symbol AS mc3_hugo_symbol, gdc.hugo_symbol AS gdc_hugo_symbol,
    gdc.reference_allele AS gdc_reference_allele, 
    m.reference_allele AS mc3_reference_allele,
    gdc.gdc_validation_status AS gdc_validation_status,
    gdc.all_effects AS gdc_all_effects, 
    m.all_effects AS mc3_all_effects,
    gdc.mc3_overlap,
    m.variant_classification AS mc3_variant_classification
FROM gdc_recoverable_unique g
LEFT JOIN mc3_protected m ON g.mc3_protected_rowid=m.rowid
LEFT JOIN gdc_shared_samples gdc ON g.gdc_rowid=gdc.rowid
;

--- MC3 recoverable unique
.once '| gzip -c > processed_data/mc3_recoverable_unique_variants.tsv.gz'
SELECT m.*,
    o.mc3_hugo_symbol, g.hugo_symbol AS gdc_hugo_symbol,
    o.reference_allele AS reference_allele,
    g.mutation_status, g.gdc_validation_status,
    g.mc3_overlap, g.somatic,
    g.dbsnp_rs AS gdc_dbsnp_rs
FROM mc3_recoverable_unique m
LEFT JOIN gdc_protected_loose_grouped g ON m.gdc_protected_rowid=g.rowid
LEFT JOIN full_overlap o ON m.overlap_rowid=o.rowid
;

-- GDC not recoverable unique
.once '| gzip -c > processed_data/gdc_not_recoverable_unique_variants.tsv.gz'
SELECT *
FROM full_overlap fo
WHERE fo.rowid NOT IN (SELECT DISTINCT overlap_rowid FROM gdc_recoverable_unique)
  AND fo.only_in_gdc = 1
;

--- MC3 not recoverable unique
.once '| gzip -c > processed_data/mc3_not_recoverable_unique_variants.tsv.gz'
SELECT *
FROM full_overlap fo
WHERE fo.rowid NOT IN (SELECT DISTINCT overlap_rowid FROM mc3_recoverable_unique)
  AND fo.only_in_mc3 = 1
;
