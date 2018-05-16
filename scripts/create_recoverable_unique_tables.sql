PRAGMA cache_size=-32000000;
PRAGMA temp_store=MEMORY;

CREATE TEMPORARY TABLE gdc_unique AS
    SELECT chromosome, start_position, end_position, tumor_sample_barcode,
           gdc_tumor_seq_allele1, tumor_seq_allele2,
           gdc_variant_type, gdc_variant_classification,
           gdc_t_depth_per_caller, gdc_t_ref_count_per_caller, gdc_t_alt_count_per_caller,
           gdc_n_depth_per_caller, gdc_n_ref_count_per_caller, gdc_n_alt_count_per_caller,
           gdc_filter, gdc_gdc_filter, gdc_callers,
           cancer_type,
           gdc_rowid, rowid AS overlap_rowid
    FROM full_overlap
    WHERE only_in_gdc = 1
;

CREATE TEMPORARY TABLE mc3_unique AS
    SELECT chromosome, start_position, end_position, tumor_sample_barcode,
           mc3_tumor_seq_allele1, tumor_seq_allele2,
           mc3_variant_type, mc3_variant_classification,
           mc3_t_depth, mc3_t_ref_count, mc3_t_alt_count,
           mc3_n_depth, mc3_n_ref_count, mc3_n_alt_count,
           mc3_filter, mc3_callers,
           cancer_type,
           mc3_rowid, rowid AS overlap_rowid
    FROM full_overlap
    WHERE only_in_mc3 = 1
;

CREATE TABLE IF NOT EXISTS gdc_recoverable_unique AS
SELECT
    gu.*,
    m.tumor_seq_allele1 AS mc3_tumor_seq_allele1,
    m.t_depth AS mc3_t_depth, m.t_ref_count AS mc3_t_ref_count, m.t_alt_count AS mc3_t_alt_count,
    m.n_depth AS mc3_n_depth, m.n_ref_count AS mc3_n_ref_count, m.n_alt_count AS mc3_n_alt_count,
    m.centers AS mc3_callers, m.ncallers AS mc3_ncallers,
    (CASE WHEN m.existing_variation='.' THEN NULL ELSE m.existing_variation END) AS mc3_existing_variation,
    m.rowid AS mc3_protected_rowid
FROM gdc_unique gu
INNER JOIN mc3_protected m
    USING (tumor_sample_barcode, chromosome, start_position, end_position, tumor_seq_allele2)
;

CREATE TABLE IF NOT EXISTS mc3_recoverable_unique AS
SELECT 
    mu.*
    g.tumor_seq_allele1 AS gdc_tumor_seq_allele1,
    g.t_depth_per_caller AS gdc_t_depth_per_caller,
    g.t_ref_count_per_caller AS gdc_t_ref_count_per_caller,
    g.t_alt_count_per_caller AS gdc_t_alt_count_per_caller,
    g.n_depth_per_caller AS gdc_n_depth_per_caller,
    g.n_ref_count_per_caller AS gdc_n_ref_count_per_caller,
    g.n_alt_count_per_caller AS gdc_n_alt_count_per_caller,
    g.filter AS gdc_filter, g.gdc_filter AS gdc_gdc_filter, g.callers AS gdc_callers,
    (CASE WHEN g.existing_variation='.' THEN NULL ELSE g.existing_variation END) AS gdc_existing_variation,
    g.rowid AS gdc_protected_rowid
FROM mc3_unique mu
INNER JOIN gdc_protected_loose_grouped g
    USING (tumor_sample_barcode, chromosome, start_position, end_position, tumor_seq_allele2)
;
