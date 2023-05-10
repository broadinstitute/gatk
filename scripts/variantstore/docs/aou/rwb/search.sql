EXPORT DATA OPTIONS (
    uri ='gs://bucket/path/to/rwb-export/rwb-export-*.tsv', format ='CSV', header = false, field_delimiter ="\t") as
    SELECT (case cast(aa.location / 1000000000000 as int64)
                when 23 then 'X'
                when 24 then 'Y'
                else cast(cast(aa.location / 1000000000000 as int64) as string) end) as chromosome,
           MOD(aa.location, 1000000000000)                                           as position,
           aa.ref,
           aa.allele,
           aa.sample_id,
           si.sample_name
    FROM `gvs-internal.gg_quickstart1.alt_allele` as aa
             join `gvs-internal.gg_quickstart1.sample_info` as si
                  on aa.sample_id = si.sample_id
             join
         (SELECT ((case SPLIT(vid, '-')[OFFSET(0)]
                       when 'X' then 23
                       when 'Y' then 24
                       else cast(SPLIT(vid, '-')[OFFSET(0)] AS int64) end) * 1000000000000 +
                  CAST(SPLIT(vid, '-')[OFFSET(1)] AS int64)) as location,
                 SPLIT(vid, '-')[OFFSET(2)]                  as ref,
                 SPLIT(vid, '-')[OFFSET(3)]                  as alt
          FROM `gvs-internal.gg_quickstart1.gg-quickstart1_vat_12`)
             as vat
         on vat.ref = aa.ref and vat.alt = aa.allele and vat.location = aa.location
    where aa.location > 20000000000000
      and aa.location < 21000000000000
    order by aa.location
