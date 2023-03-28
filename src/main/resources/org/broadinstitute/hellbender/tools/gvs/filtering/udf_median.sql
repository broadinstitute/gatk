--    This is a local implementation of median, as defined in the open source bigquery udf library.  It was previously
--    referenced at "`bqutil`.fn.median"
--    However, bqutil is not installed in locations that are not "US" so we will define it ourselves.
--    https://github.com/GoogleCloudPlatform/bigquery-utils/blob/master/udfs/community/median.sqlx

CREATE TEMPORARY FUNCTION median(arr ANY TYPE) AS ((
  SELECT IF (
    MOD(ARRAY_LENGTH(arr), 2) = 0,
    (arr[OFFSET(DIV(ARRAY_LENGTH(arr), 2) - 1)] + arr[OFFSET(DIV(ARRAY_LENGTH(arr), 2))]) / 2,
    arr[OFFSET(DIV(ARRAY_LENGTH(arr), 2))]
  )
  FROM (SELECT ARRAY_AGG(x ORDER BY x) AS arr FROM UNNEST(arr) AS x)
));
