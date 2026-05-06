-- This query identifies loci at which VIDs present in the R2 dataset but not in the classic dataset.
-- It extracts CHROM-POS from the VID and counts the occurrences of each, grouping and ordering by descending count.
SELECT regexp_extract(vid, '^([^-]+-[0-9]+)'), COUNT(*)
FROM
    (
        SELECT DISTINCT r2.vid
        FROM foxtrot.foxtrot_v4_2025_07_29_vat_v9_r2 r2
                 LEFT OUTER JOIN
             foxtrot.foxtrot_v4_2025_07_29_vat classic
             ON r2.vid = classic.vid
        WHERE classic.vid IS NULL
    )
GROUP BY 1
ORDER BY 2 DESC
