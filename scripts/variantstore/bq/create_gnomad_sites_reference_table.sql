CREATE TABLE `spec-ops-aou.gvs_public_reference_data.gnomad_v3_sites` AS
SELECT chrom * 1000000000000 + start_position as location
FROM (
    SELECT CASE WHEN reference_name = "chrX" THEN 23 ELSE CASE WHEN reference_name = "chrY" THEN 24 ELSE CAST(REPLACE(reference_name,"chr","") AS INT64) END END chrom, start_position
    FROM (
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr1`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr2`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr3`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr4`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr5`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr6`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr7`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr8`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr9`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr10`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr11`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr12`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr13`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr14`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr15`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr16`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr17`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr18`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr19`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr20`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr21`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chr22`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chrX`
        UNION ALL
        SELECT * FROM `bigquery-public-data.gnomAD.v3_genomes__chrY`        
    )
)  
