PRAGMA cache_size=-4192000;
PRAGMA temp_store=MEMORY;
PRAGMA journal_mode=OFF;

-- Subset samples to be shared from MC3 and GDC
-- MC3
DROP TABLE IF EXISTS mc3_shared_samples;

CREATE TABLE IF NOT EXISTS mc3_shared_samples AS
SELECT * FROM mc3
WHERE tumor_sample_barcode IN (
    SELECT DISTINCT tumor_sample_barcode FROM gdc_grouped_callers
);
CREATE INDEX ix_mc3_shared_samples_tumor_barcode ON mc3_shared_samples (tumor_sample_barcode);


-- GDC
DROP TABLE IF EXISTS gdc_shared_samples;

CREATE TABLE IF NOT EXISTS gdc_shared_samples AS
SELECT * FROM gdc_grouped_callers
WHERE tumor_sample_barcode IN (
    SELECT DISTINCT tumor_sample_barcode FROM mc3_shared_samples
    INTERSECT
    SELECT DISTINCT tumor_sample_barcode FROM gdc_grouped_callers
);
CREATE INDEX ix_gdc_shared_samples_tumor_barcode ON gdc_shared_samples (tumor_sample_barcode);

