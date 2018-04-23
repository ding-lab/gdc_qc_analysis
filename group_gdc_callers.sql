DROP TABLE IF EXISTS gdc_grouped_callers;
CREATE TABLE IF NOT EXISTS gdc_grouped_callers AS
SELECT
    hugo_symbol, entrez_gene_id, center, ncbi_build,
    chromosome, start_position, end_position, strand,
    variant_classification, variant_type,
    reference_allele, tumor_seq_allele1, tumor_seq_allele2,
    dbsnp_rs, dbsnp_val_status,
    tumor_sample_barcode, matched_norm_sample_barcode,
    match_norm_seq_allele1, match_norm_seq_allele2,
    tumor_validation_allele1, tumor_validation_allele2,
    match_norm_validation_allele1, match_norm_validation_allele2, verification_status, validation_status, mutation_status, sequencing_phase,
    sequence_source, validation_method, score, bam_file, sequencer,
    tumor_sample_uuid, matched_norm_sample_uuid, hgvsc, hgvsp, hgvsp_short,
    transcript_id, exon_number,
    all_effects, allele,
    group_concat(t_depth, ',') AS t_depth_per_caller, group_concat(t_ref_count, ',') AS t_ref_count_per_caller, group_concat(t_alt_count, ',') AS t_alt_count_per_caller,
    group_concat(n_depth, ',') AS n_depth_per_caller, group_concat(n_ref_count, ',') AS n_ref_count_per_caller, group_concat(n_alt_count, ',') AS n_alt_count_per_caller,
    gene, feature, feature_type, one_consequence, consequence,
    cdna_position, cds_position, protein_position,
    amino_acids, codons, existing_variation, allele_num, distance,
    transcript_strand, symbol, symbol_source, hgnc_id, biotype,
    canonical, ccds, ensp, swissprot, trembl, uniparc, refseq,
    sift, polyphen, exon, intron, domains,
    gmaf, afr_maf, amr_maf, asn_maf, eas_maf, eur_maf, sas_maf, aa_maf, ea_maf,
    clin_sig, somatic, pubmed, motif_name, motif_pos, high_inf_pos,
    motif_score_change, impact, pick, variant_class, tsl, hgvs_offset, pheno,
    minimised, exac_af, exac_af_adj, exac_af_afr, exac_af_amr, exac_af_eas, exac_af_fin, exac_af_nfe, exac_af_oth, exac_af_sas,
    gene_pheno, 'filter', context, src_vcf_id, tumor_bam_uuid, normal_bam_uuid, case_id, gdc_filter, cosmic, mc3_overlap, gdc_validation_status,
	cancer_type, group_concat(caller, '|') AS callers
FROM gdc
GROUP BY
    chromosome, start_position, end_position, strand,
    tumor_seq_allele1, tumor_seq_allele2,
    tumor_sample_barcode, matched_norm_sample_barcode
ORDER BY tumor_sample_barcode, chromosome, start_position;
