library(tidyverse)

full_tbl <- read_tsv(
    '../GDC_QC_data/overlap/fullouterjoin.tsv',
    col_types = cols(
        .default = col_character(),
        entrez_gene_id = col_integer(),
        start_position = col_integer(),
        end_position = col_integer(),
        t_depth = col_integer(),
        t_ref_count = col_integer(),
        t_alt_count = col_integer(),
        n_depth = col_integer(),
        n_ref_count = col_integer(),
        n_alt_count = col_integer(),
        strand_vep = col_integer(),
        ncallers = col_integer(),
        `entrez_gene_id:1` = col_integer(),
        `start_position:1` = col_integer(),
        `end_position:1` = col_integer(),
        t_depth_per_caller = col_number(),
        t_ref_count_per_caller = col_number(),
        t_alt_count_per_caller = col_number(),
        n_depth_per_caller = col_number(),
        `allele_num:1` = col_integer(),
        `distance:1` = col_integer(),
        transcript_strand = col_integer(),
        `gmaf:1` = col_double(),
        `afr_maf:1` = col_double(),
        `amr_maf:1` = col_double(),
        `eas_maf:1` = col_double(),
        `eur_maf:1` = col_double(),
        `sas_maf:1` = col_double(),
        `aa_maf:1` = col_double(),
        `ea_maf:1` = col_double(),
        `pick:1` = col_integer(),
        `tsl:1` = col_integer(),
        `hgvs_offset:1` = col_integer(),
        `minimised:1` = col_integer(),
        `exac_af:1` = col_double(),
        exac_af_adj = col_double(),
        `exac_af_afr:1` = col_double(),
        `exac_af_amr:1` = col_double(),
        `exac_af_eas:1` = col_double(),
        `exac_af_fin:1` = col_double(),
        `exac_af_nfe:1` = col_double(),
        `exac_af_oth:1` = col_double(),
        `exac_af_sas:1` = col_double()
    ))

# Fill in tumor sample barcode if it's NA
full_tbl$tumor_sample_barcode[is.na(full_tbl$tumor_sample_barcode)] <- full_tbl$`tumor_sample_barcode:1`[is.na(full_tbl$tumor_sample_barcode)]

sample_cancer_type_tbl <- full_tbl %>% 
    select(tumor_sample_barcode, cancer_type) %>% 
    filter(!is.na(cancer_type)) %>% 
    distinct()

grp_sample <- group_by(full_tbl, tumor_sample_barcode)
per_sample_overlap <- grp_sample %>%
    summarize(
        n_shared = sum(!is.na(centers) & !is.na(callers)),
        n_mc3_unique = sum(!is.na(centers) & is.na(callers)),
        n_gdc_unique = sum(is.na(centers) & !is.na(callers)),
        n_total= n()
    ) %>%
    mutate(
        mc3_shared_percent = n_shared / (n_mc3_unique + n_shared),
        gdc_shared_percent = n_shared / (n_gdc_unique + n_shared)
    ) %>% 
    left_join(sample_cancer_type_tbl, by = "tumor_sample_barcode") %>%
    arrange(cancer_type)

write_tsv(per_sample_overlap, '~/Box Sync/Ding_Lab/Projects_Current/GDC-QC/Results/per_sample_overlap.tsv')

ggplot(per_sample_overlap, aes(x=gdc_shared_percent, y=mc3_shared_percent, size=n_total)) +
    geom_point(alpha=0.65) + 
    scale_x_continuous(name="GDC shared calls", labels = scales::percent) +
    scale_y_continuous(name="MC3 shared calls", labels = scales::percent) + 
    scale_color_discrete(name="Cancer type") + 
    scale_size_area(name="Total calls", max_size = 5, breaks = c(10, 100, 1000, 10000), trans='sqrt') + 
    facet_wrap(~ cancer_type, nrow = 2) +
    theme_bw() + 
    labs(title = "Somatic variant call overlap per sample")

# BRCA bottom right two samples
per_sample_overlap %>%
    filter(cancer_type == 'BRCA' & gdc_shared_percent > 0.6 & mc3_shared_percent < 0.12)

# BRCA top left one sample
per_sample_overlap %>%
    filter(cancer_type == 'BRCA' & gdc_shared_percent < 0.12 & mc3_shared_percent > 0.8)

# COAD top left two samples
per_sample_overlap %>%
    filter(cancer_type == 'COAD' & gdc_shared_percent < 0.15 & mc3_shared_percent > 0.8)

# COAD bottom right two samples
per_sample_overlap %>%
    filter(cancer_type == 'COAD' & gdc_shared_percent > 0.75 & mc3_shared_percent < 0.35)

# OV middle left biggest sample
per_sample_overlap %>%
    filter(cancer_type == 'OV' & gdc_shared_percent < 0.14 & mc3_shared_percent > 0.5) %>%
    arrange(-n_total) %>%
    head(1)

# OV bottom right two samples
per_sample_overlap %>%
    filter(cancer_type == 'OV' & gdc_shared_percent > 0.48 & mc3_shared_percent < 0.10)

# BRCA middle right group of samples
brca_samples <- per_sample_overlap %>%
    filter(cancer_type == 'BRCA' & mc3_shared_percent > 0.37 & mc3_shared_percent < 0.55 & 
               gdc_shared_percent > 0.5 & gdc_shared_percent < 1) %>%
    arrange(mc3_shared_percent)

# COAD middle right group of samples
coad_samples <- per_sample_overlap %>%
    filter(cancer_type == 'COAD' & mc3_shared_percent > 0.4 & mc3_shared_percent < 0.6 &
               gdc_shared_percent > 0.7 & gdc_shared_percent < 1) %>%
    arrange(mc3_shared_percent)

brca_tbl <- full_tbl %>% 
    filter(tumor_sample_barcode %in% brca_samples$tumor_sample_barcode) 
coad_tbl <- full_tbl %>% 
    filter(tumor_sample_barcode %in% coad_samples$tumor_sample_barcode)

# Check sample center
# Center code: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes
table(str_sub(brca_samples$tumor_sample_barcode, 27))  # 09 is WUGSC
table(str_sub(coad_samples$tumor_sample_barcode, 27))  # 10 is BCM

# Check plate ID
table(str_sub(brca_samples$tumor_sample_barcode, 22, 25)) %>% as_tibble %>% arrange(-n)
table(str_sub(coad_samples$tumor_sample_barcode, 22, 25)) %>% as_tibble %>% arrange(-n)

# Check filters
brca_mc3_unique <- brca_tbl %>% 
    filter(!is.na(centers) & is.na(callers)) %>%
    group_by(`filter`, `filter:1`) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(-n)
brca_mc3_unique %>% 
    filter(str_detect(filter, 'nonpreferredpair')) %>% 
    summarise(total = sum(n))
brca_samples %>% summarise(total_mc3_unique = sum(n_mc3_unique))


# Check filters
coad_mc3_unique <- coad_tbl %>% 
    filter(!is.na(centers) & is.na(callers)) %>%
    group_by(`filter`, `filter:1`) %>%
    summarise(n = n()) %>%
    arrange(-n) %>%
    ungroup()

coad_mc3_unique %>% 
    filter(str_detect(filter, 'nonpreferredpair')) %>% 
    summarise(total = sum(n))
coad_samples %>% summarise(total_mc3_unique = sum(n_mc3_unique))

# BRCA top left
brca_samples <- per_sample_overlap %>%
    filter(cancer_type == 'BRCA' & mc3_shared_percent > 0.8 & gdc_shared_percent < 0.1)
brca_tbl <- full_tbl %>% 
    filter(tumor_sample_barcode %in% brca_samples$tumor_sample_barcode)
brca_tbl %>% 
    filter(is.na(centers) & !is.na(callers)) %>%
    group_by(gdc_validation_status) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(-n)

# COAD top left
coad_samples <- per_sample_overlap %>%
    filter(cancer_type == 'BRCA' & mc3_shared_percent > 0.8 & gdc_shared_percent < 0.25)
coad_tbl <- full_tbl %>% 
    filter(tumor_sample_barcode %in% coad_samples$tumor_sample_barcode)
coad_tbl %>% 
    filter(is.na(centers) & !is.na(callers)) %>%
    group_by(callers) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(-n)

# OV bottom right
ov_samples <- per_sample_overlap %>%
    filter(cancer_type == 'OV' & mc3_shared_percent < .1 & gdc_shared_percent > 0.25)
ov_tbl <- full_tbl %>% 
    filter(tumor_sample_barcode %in% ov_samples$tumor_sample_barcode)
ov_tbl %>% 
    filter(!is.na(centers) & is.na(callers)) %>%
    group_by(filter) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(-n)


# Group GDC unique calls
gdc_calls_tbl <- full_tbl %>% 
    filter(!is.na(callers)) 

gdc_dup_tbl <- gdc_calls_tbl %>%
    select(
        tumor_sample_barcode, `chromosome:1`, `start_position:1`, `end_position:1`, 
        `reference_allele:1`,`tumor_seq_allele1:1`, `tumor_seq_allele2:1`,
        callers,
        `filter:1`, gdc_filter, gdc_validation_status, `validation_method:1`, `validation_status:1`,
        `filter`, `hugo_symbol:1`
    ) %>%
    distinct(
        tumor_sample_barcode, `chromosome:1`, `start_position:1`, `end_position:1`, 
        `tumor_seq_allele1:1`, `tumor_seq_allele2:1`, 
        .keep_all = TRUE
    ) %>%
    group_by(`chromosome:1`, `start_position:1`, `end_position:1`, tumor_sample_barcode) %>%
    filter(n() > 1)

gdc_dup_tbl %>% 
    arrange(`chromosome:1`, `start_position:1`, `end_position:1`, tumor_sample_barcode) %>%
    View


gdc_dup_tbl %>%
    group_by(`chromosome:1`, `start_position:1`, `end_position:1`) %>%
    filter(n() >= 10) %>%
    arrange(`chromosome:1`, `start_position:1`, `end_position:1`, tumor_sample_barcode) %>%
    View

