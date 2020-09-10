WITH oper4 AS (
  SELECT 41200 FieldID, eid, oper4 code, 
    CASE 
      WHEN h.admidate IS NOT NULL THEN h.admidate
      WHEN h.admidate IS NULL AND h.opdate IS NOT NULL THEN h.opdate 
      ELSE h.epistart
    END vdate
  FROM `broad-ml4cvd.ukbb7089_201904.hesin` h
  WHERE oper4 IS NOT NULL
), diag_icd10 AS (
  SELECT 41202 FieldID, eid, diag_icd10 code,
    CASE 
      WHEN h.admidate IS NOT NULL THEN h.admidate
      WHEN h.admidate IS NULL AND h.opdate IS NOT NULL THEN h.opdate 
      ELSE h.epistart
    END vdate
  FROM `broad-ml4cvd.ukbb7089_201904.hesin` h
  WHERE diag_icd10 IS NOT NULL
), diag_icd9 AS (
  SELECT 41203 FieldID, eid, diag_icd9 code,
    CASE 
      WHEN h.admidate IS NOT NULL THEN h.admidate
      WHEN h.admidate IS NULL AND h.opdate IS NOT NULL THEN h.opdate 
      ELSE h.epistart
    END vdate
  FROM `broad-ml4cvd.ukbb7089_201904.hesin` h
  WHERE diag_icd9 IS NOT NULL
), oper4secondary AS (
  SELECT 41210 FieldID, h.eid, sec.oper4 code, 
    CASE 
      WHEN h.admidate IS NOT NULL THEN h.admidate
      WHEN h.admidate IS NULL AND h.opdate IS NOT NULL THEN h.opdate 
      ELSE h.epistart
    END vdate
  FROM `broad-ml4cvd.ukbb7089_201904.hesin_oper` sec
  LEFT JOIN `broad-ml4cvd.ukbb7089_201904.hesin` h ON sec.eid=h.eid AND sec.record_id=h.record_id
  WHERE TRUE
    AND sec.oper4 IS NOT NULL
), diag_icd10_secondary AS (
  SELECT 41204 FieldID, h.eid, sec.diag_icd10 code, 
    CASE 
      WHEN h.admidate IS NOT NULL THEN h.admidate
      WHEN h.admidate IS NULL AND h.opdate IS NOT NULL THEN h.opdate 
      ELSE h.epistart
    END vdate
  FROM `broad-ml4cvd.ukbb7089_201904.hesin_diag10` sec
  LEFT JOIN `broad-ml4cvd.ukbb7089_201904.hesin` h ON sec.eid=h.eid AND sec.record_id=h.record_id
  WHERE TRUE
    AND sec.diag_icd10 IS NOT NULL
), diag_icd9_secondary AS (
  SELECT 41205 FieldID, h.eid, sec.diag_icd9 code, 
    CASE 
      WHEN h.admidate IS NOT NULL THEN h.admidate
      WHEN h.admidate IS NULL AND h.opdate IS NOT NULL THEN h.opdate 
      ELSE h.epistart
    END vdate
  FROM `broad-ml4cvd.ukbb7089_201904.hesin_diag9` sec
  LEFT JOIN `broad-ml4cvd.ukbb7089_201904.hesin` h ON sec.eid=h.eid AND sec.record_id=h.record_id
  WHERE TRUE
    AND sec.diag_icd9 IS NOT NULL
)

SELECT 
  diagnostics.eid sample_id, diagnostics.FieldID, diagnostics.code value, 
  CASE 
    WHEN MIN(PARSE_DATE("%E4Y-%m-%d", vdate)) IS NULL THEN MIN(PARSE_DATE("%E4Y-%m-%d", p.value))
    ELSE MIN(PARSE_DATE("%E4Y-%m-%d", vdate))
  END first_date
FROM (
  SELECT * FROM oper4
  UNION DISTINCT
  SELECT * FROM diag_icd10
  UNION DISTINCT
  SELECT * FROM diag_icd9
  UNION DISTINCT
  SELECT * FROM oper4secondary
  UNION DISTINCT
  SELECT * FROM diag_icd10_secondary
  UNION DISTINCT
  SELECT * FROM diag_icd9_secondary
) diagnostics
JOIN `broad-ml4cvd.ukbb7089_201904.phenotype` p ON p.sample_id = diagnostics.eid AND p.array_idx=0 AND p.instance=0 AND p.FieldID=53
GROUP BY diagnostics.eid, diagnostics.FieldID, diagnostics.code
ORDER BY first_date ASC
;