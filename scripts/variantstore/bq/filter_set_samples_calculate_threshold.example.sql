WITH fss AS (
  SELECT *, 
         (ins_count / del_count) as ins_del_ratio, 
         (ti_count / tv_count) as ti_tv_ratio, 
         (snp_het_count / snp_homvar_count) snp_het_homvar_ratio, 
         (indel_het_count / indel_homvar_count) as indel_het_homvar_ratio
  FROM `aou-genomics-curation-prod.alpha1_1000.filter_set_samples`
  WHERE filter_set_name = 'kc-testing-1'),
medians AS ( 
    SELECT 
        `bqutil`.fn.median(ARRAY_AGG(del_count IGNORE NULLS)) as m_del_count,
        `bqutil`.fn.median(ARRAY_AGG(ins_count IGNORE NULLS)) as m_ins_count,
        `bqutil`.fn.median(ARRAY_AGG(snp_count IGNORE NULLS)) as m_snp_count,
        `bqutil`.fn.median(ARRAY_AGG(singleton IGNORE NULLS)) as m_singleton,
        `bqutil`.fn.median(ARRAY_AGG(ins_del_ratio IGNORE NULLS)) as m_ins_del_ratio,
        `bqutil`.fn.median(ARRAY_AGG(ti_tv_ratio IGNORE NULLS)) as m_ti_tv_ratio,
        `bqutil`.fn.median(ARRAY_AGG(snp_het_homvar_ratio IGNORE NULLS)) as m_snp_het_homvar_ratio,
        `bqutil`.fn.median(ARRAY_AGG(indel_het_homvar_ratio IGNORE NULLS)) as m_indel_het_homvar_ratio        
    FROM fss),
mads AS (
    SELECT 
        `bqutil`.fn.median(ARRAY_AGG(ABS(del_count - m_del_count) IGNORE NULLS)) as mad_del_count,
        `bqutil`.fn.median(ARRAY_AGG(ABS(ins_count - m_ins_count) IGNORE NULLS)) as mad_ins_count,
        `bqutil`.fn.median(ARRAY_AGG(ABS(snp_count - m_snp_count) IGNORE NULLS)) as mad_snp_count,
        `bqutil`.fn.median(ARRAY_AGG(ABS(singleton - m_singleton) IGNORE NULLS)) as mad_singleton,
        `bqutil`.fn.median(ARRAY_AGG(ABS(ins_del_ratio - m_ins_del_ratio) IGNORE NULLS)) as mad_ins_del_ratio,
        `bqutil`.fn.median(ARRAY_AGG(ABS(ti_tv_ratio - m_ti_tv_ratio) IGNORE NULLS)) as mad_ti_tv_ratio,
        `bqutil`.fn.median(ARRAY_AGG(ABS(snp_het_homvar_ratio - m_snp_het_homvar_ratio) IGNORE NULLS)) as mad_snp_het_homvar_ratio,
        `bqutil`.fn.median(ARRAY_AGG(ABS(indel_het_homvar_ratio - m_indel_het_homvar_ratio) IGNORE NULLS)) as mad_indel_het_homvar_ratio
    FROM fss
    CROSS JOIN medians 
    WHERE filter_set_name = 'kc-testing-1')

SELECT 
    fss.sample_id,
    si.sample_name,
    del_count, m_del_count, mad_del_count,
    CASE WHEN del_count BETWEEN m_del_count - 4*mad_del_count AND m_del_count + 4*mad_del_count THEN true ELSE false END pass_del_count,
    
    ins_count, m_ins_count, mad_ins_count,
    CASE WHEN ins_count BETWEEN m_ins_count - 4*mad_ins_count AND m_ins_count + 4*mad_ins_count THEN true ELSE false END pass_ins_count,

    snp_count, m_snp_count, mad_snp_count,
    CASE WHEN snp_count BETWEEN m_snp_count - 4*mad_snp_count AND m_snp_count + 4*mad_snp_count THEN true ELSE false END pass_snp_count,

    singleton, m_singleton, mad_singleton,
    CASE WHEN singleton BETWEEN m_singleton - 8*mad_singleton AND m_singleton + 8*mad_singleton THEN true ELSE false END pass_singleton,
    
    ins_del_ratio, m_ins_del_ratio, mad_ins_del_ratio,
    CASE WHEN ins_del_ratio BETWEEN m_ins_del_ratio - 4*mad_ins_del_ratio AND m_ins_del_ratio + 4*mad_ins_del_ratio THEN true ELSE false END pass_ins_del_ratio,

    ti_tv_ratio, m_ti_tv_ratio, mad_ti_tv_ratio,
    CASE WHEN ti_tv_ratio BETWEEN m_ti_tv_ratio - 4*mad_ti_tv_ratio AND m_ti_tv_ratio + 4*mad_ti_tv_ratio THEN true ELSE false END pass_ti_tv_ratio,

    snp_het_homvar_ratio, m_snp_het_homvar_ratio, mad_snp_het_homvar_ratio,
    CASE WHEN snp_het_homvar_ratio BETWEEN m_snp_het_homvar_ratio - 4*mad_snp_het_homvar_ratio AND m_snp_het_homvar_ratio + 4*mad_snp_het_homvar_ratio THEN true ELSE false END pass_snp_het_homvar_ratio,

    indel_het_homvar_ratio, m_indel_het_homvar_ratio, mad_indel_het_homvar_ratio,
    CASE WHEN indel_het_homvar_ratio BETWEEN m_indel_het_homvar_ratio - 4*mad_indel_het_homvar_ratio AND m_indel_het_homvar_ratio + 4*mad_indel_het_homvar_ratio THEN true ELSE false END pass_indel_het_homvar_ratio,    
FROM fss
JOIN `aou-genomics-curation-prod.alpha1_1000.sample_info` si ON (fss.sample_id = si.sample_id)
CROSS JOIN medians 
CROSS JOIN mads
order by 1