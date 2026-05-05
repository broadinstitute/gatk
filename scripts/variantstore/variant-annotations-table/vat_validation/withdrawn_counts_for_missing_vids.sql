-- Query to identify non-control samples associated with VIDs that were in the classic VAT
-- but which are not in the r2 VAT. Reports the withdrawn date and the number of samples
-- correlating with that date.
SELECT withdrawn, COUNT(*) AS count
FROM `foxtrot.sample_info` si
WHERE
  si.is_control IS FALSE
  AND CAST(si.sample_name AS INT64)
    IN (
      SELECT DISTINCT
        person_id
      FROM
        `foxtrot.vid_to_participant_mapping_2025_10_28`
          a
      CROSS JOIN
        UNNEST(a.person_ids) AS person_id
      WHERE
        NOT EXISTS(
          SELECT 1
          FROM
            `foxtrot.foxtrot_v4_2025_07_29_vat_v9_r2`
              b
          WHERE a.vid = b.vid
        )
      ORDER BY person_id
    )
GROUP BY withdrawn
ORDER BY withdrawn
