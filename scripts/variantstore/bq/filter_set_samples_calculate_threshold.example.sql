WITH medians AS ( 
    SELECT 
        `bqutil`.fn.median(ARRAY_AGG(del_count IGNORE NULLS)) as m_del_count,
        `bqutil`.fn.median(ARRAY_AGG(ins_count IGNORE NULLS)) as m_ins_count,
        `bqutil`.fn.median(ARRAY_AGG(snp_count IGNORE NULLS)) as m_snp_count,
        `bqutil`.fn.median(ARRAY_AGG(singleton IGNORE NULLS)) as m_singleton
    FROM `aou-genomics-curation-prod.alpha1_1000.filter_set_samples` 
    WHERE filter_set_name = 'kc-testing-1'),
mads AS (
    SELECT 
        `bqutil`.fn.median(ARRAY_AGG(ABS(del_count - m_del_count) IGNORE NULLS)) as mad_del_count,
        `bqutil`.fn.median(ARRAY_AGG(ABS(ins_count - m_ins_count) IGNORE NULLS)) as mad_ins_count,
        `bqutil`.fn.median(ARRAY_AGG(ABS(snp_count - m_snp_count) IGNORE NULLS)) as mad_snp_count,
        `bqutil`.fn.median(ARRAY_AGG(ABS(singleton - m_singleton) IGNORE NULLS)) as mad_singleton
    FROM `aou-genomics-curation-prod.alpha1_1000.filter_set_samples`
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
    CASE WHEN singleton BETWEEN m_singleton - 8*mad_singleton AND m_singleton + 8*mad_singleton THEN true ELSE false END pass_singleton
FROM `aou-genomics-curation-prod.alpha1_1000.filter_set_samples` fss
JOIN `aou-genomics-curation-prod.alpha1_1000.sample_info` si ON (fss.sample_id = si.sample_id)
CROSS JOIN medians 
CROSS JOIN mads
WHERE filter_set_name = 'kc-testing-1'
order by 1