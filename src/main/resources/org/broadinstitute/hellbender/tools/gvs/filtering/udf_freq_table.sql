CREATE TEMPORARY FUNCTION freq_table(arr ANY TYPE) AS (
( -- need these extra parentheses to use this query as an expression
SELECT ARRAY_AGG(entry)
FROM (
SELECT STRUCT(x as value, COUNT(1) as freq) as entry
FROM UNNEST(arr) AS x
GROUP BY x
ORDER BY x
)
)
);
