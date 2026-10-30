-- Verification harness for the participant mapping table scrub (VDS/person-mapping discrepancy).
-- Pure CTE fixture: creates nothing, reads no real table. Safe to run in any project.
--   bq query --nouse_legacy_sql --project_id=<any> < verify_participant_mapping_scrub.sql
-- Every row of the result must read 'pass'.

WITH
-- ---------- fixture: sample_info ----------
sample_info AS (
  SELECT * FROM UNNEST([
    STRUCT('1001'    AS sample_name, 1 AS sample_id, false AS is_control, CAST(NULL AS TIMESTAMP)          AS withdrawn),
           ('1002',                  2,              false,               TIMESTAMP '2025-07-29 00:00:00+00'),
           ('1003',                  3,              false,               TIMESTAMP '2025-07-29 00:00:00+00'),
           ('1004',                  4,              true,                CAST(NULL AS TIMESTAMP)),
           ('NA12878',               5,              true,                CAST(NULL AS TIMESTAMP)),
           ('1006',                  6,              false,               TIMESTAMP '2025-07-29 00:00:00+00'),
           ('1006',                  7,              false,               CAST(NULL AS TIMESTAMP)),
           ('SOMENAME',               8,              false,               CAST(NULL AS TIMESTAMP))
  ])
),
-- ---------- fixture: participant mapping table as delivered ----------
participant_mapping AS (
  SELECT * FROM UNNEST([
    STRUCT('v_allgood'  AS vid, [1001]                   AS person_ids),
           ('v_mixed',           [1001, 1002, 1003]),
           ('v_allbad',          [1002, 1003]),
           ('v_control',         [1001, 1004]),
           ('v_reingest',        [1006]),
           ('v_order',           [1003, 1001, 1002, 1006]),
           ('v_dup',             [1001, 1001, 1002]),
           ('v_empty',           CAST([] AS ARRAY<INT64>))
  ])
),

-- ================= LOGIC UNDER TEST =================
-- A person id is good iff at least one sample_info row bearing that sample_name is
-- neither withdrawn nor a control. Stated this way, a person re-ingested under a new
-- sample_id survives even though an older withdrawn row for the same name exists.
good_person_ids AS (
  SELECT DISTINCT SAFE_CAST(sample_name AS INT64) AS person_id
  FROM sample_info
  WHERE withdrawn IS NULL
    AND is_control = false
    AND SAFE_CAST(sample_name AS INT64) IS NOT NULL
),
exploded AS (
  SELECT m.vid AS vid, p AS person_id, off
  FROM participant_mapping AS m, UNNEST(m.person_ids) AS p WITH OFFSET AS off
),
scrubbed AS (
  -- INNER JOIN drops bad person ids; the GROUP BY then drops any VID left with none.
  -- ORDER BY off preserves each surviving array's original relative order.
  SELECT e.vid AS vid, ARRAY_AGG(e.person_id ORDER BY e.off) AS person_ids
  FROM exploded AS e
  JOIN good_person_ids AS g ON e.person_id = g.person_id
  GROUP BY e.vid
),
-- ====================================================

expected AS (
  SELECT * FROM UNNEST([
    STRUCT('v_allgood'  AS vid, [1001]       AS person_ids),
           ('v_mixed',           [1001]),
           ('v_control',         [1001]),
           ('v_reingest',        [1006]),
           ('v_order',           [1001, 1006]),
           ('v_dup',             [1001, 1001])
  ])
),
-- VIDs whose every carrier was withdrawn or control must disappear from the table.
expected_dropped AS (
  SELECT dropped_vid FROM UNNEST(['v_allbad', 'v_empty']) AS dropped_vid
),
checks AS (
  SELECT
    COALESCE(s.vid, x.vid) AS vid,
    TO_JSON_STRING(x.person_ids) AS expected_person_ids,
    TO_JSON_STRING(s.person_ids) AS actual_person_ids,
    CASE
      WHEN s.vid IS NULL THEN 'FAIL: row missing'
      WHEN x.vid IS NULL THEN 'FAIL: unexpected row'
      WHEN TO_JSON_STRING(s.person_ids) = TO_JSON_STRING(x.person_ids) THEN 'pass'
      ELSE 'FAIL: array mismatch'
    END AS result
  FROM scrubbed AS s
  FULL OUTER JOIN expected AS x ON s.vid = x.vid

  UNION ALL

  -- Asserted separately: the FULL OUTER JOIN above cannot express these, because a VID
  -- absent from both sides yields no row at all, making "correctly dropped"
  -- indistinguishable from "never tested".
  SELECT
    d.dropped_vid AS vid,
    '<row absent>' AS expected_person_ids,
    IFNULL((SELECT TO_JSON_STRING(s2.person_ids) FROM scrubbed AS s2 WHERE s2.vid = d.dropped_vid),
           '<row absent>') AS actual_person_ids,
    IF(EXISTS(SELECT 1 FROM scrubbed AS s2 WHERE s2.vid = d.dropped_vid),
       'FAIL: row should have been dropped', 'pass') AS result
  FROM expected_dropped AS d
)
SELECT * FROM checks ORDER BY vid
