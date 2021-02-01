CREATE TEMPORARY FUNCTION minimize(ref STRING, allele STRING)
RETURNS STRING 
  LANGUAGE js AS """
  let done = false
    while (!done && ref.length !== 1) {
        if (ref.slice(-1) === allele.slice(-1)) {
            ref = ref.slice(0, -1)
	    allele = allele.slice(0,-1)
        } else {
            done = true
        }
    }
    return ref+','+allele
    """;
    
CREATE OR REPLACE TABLE
  `spec-ops-aou.ccdg_7596_benchmark.alt_allele`
PARTITION BY
   RANGE_BUCKET(location, GENERATE_ARRAY(0, 25000000000000, 1000000000000))  
CLUSTER BY location, sample_id
AS
  -- data will be in 3 positions: 0, 1, 2
  -- we need the position 1 data from GTs that are 0/1 and 1/2
WITH
  position1 as (select * from `spec-ops-aou.ccdg_7596_benchmark.vet_001` where call_GT in ("0/1", "1/1", "0|1", "1|1")),
  -- we only need position 2 data from GTs that are 0/2
  -- QUESTION what about 0/2?
  position2 as (select * from `spec-ops-aou.ccdg_7596_benchmark.vet_001` where call_GT IN ("1/2", "1|2"))


-- we don't really need to save data for spanning deletions (i.e. allele = '*')
select location, sample_id, 
SPLIT(minimize(ref, SPLIT(alt,",")[OFFSET(0)]))[OFFSET(0)] as ref,
SPLIT(minimize(ref, SPLIT(alt,",")[OFFSET(0)]))[OFFSET(1)] as allele,
1 as allele_pos, call_GT,
as_raw_mq,
cast(SPLIT(as_raw_mq,"|")[OFFSET(1)] as int64) raw_mq,
as_raw_mqranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_mqranksum,",")[SAFE_OFFSET(0)] as float64) * 10.0 as int64) as raw_mqranksum_x_10, # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(0) for the ALT.  We also have some NaNs (safe_cast)
as_qualapprox,
cast(SPLIT(as_qualapprox,"|")[OFFSET(0)] as int64) as qual, # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(0) for the ALT
as_raw_readposranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_readposranksum,",")[SAFE_OFFSET(0)] as float64) * 10.0 as int64) as raw_readposranksum_x_10, # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(0) for the ALT.    We also have some NaNs 
#NULL as raw_readposranksum_x_10,
as_sb_table,
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(0)],",")[OFFSET(0)] as int64) as sb_ref_plus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(0)],",")[OFFSET(1)] as int64) as sb_ref_minus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(1)],",")[OFFSET(0)] as int64) as sb_alt_plus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(1)],",")[OFFSET(1)] as int64) as sb_alt_minus, 
call_AD,
cast(SPLIT(call_AD,",")[OFFSET(0)] as int64) as ref_ad, 
cast(SPLIT(call_AD,",")[OFFSET(1)] as int64) as ad
from position1

union all

select location, sample_id, 
SPLIT(minimize(ref, SPLIT(alt,",")[OFFSET(0)]))[OFFSET(0)] as ref,
SPLIT(minimize(ref, SPLIT(alt,",")[OFFSET(0)]))[OFFSET(1)] as allele,
1 as allele_pos, call_GT,
as_raw_mq,
cast(SPLIT(as_raw_mq,"|")[OFFSET(1)] as int64) raw_mq,
as_raw_mqranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_mqranksum,",")[SAFE_OFFSET(0)] as float64) * 10.0 as int64) as raw_mqranksum_x_10,  # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(0) for the ALT.  We also have some NaNs 
as_qualapprox,
cast(SPLIT(as_qualapprox,"|")[OFFSET(0)] as int64) as qual,  # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(0) for the ALT.  We also have some NaNs 
as_raw_readposranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_readposranksum,",")[SAFE_OFFSET(0)] as float64) * 10.0 as int64) as raw_readposranksum_x_10,  # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(0) for the ALT.  We also have some NaNs 
as_sb_table,
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(0)],",")[OFFSET(0)] as int64) as sb_ref_plus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(0)],",")[OFFSET(1)] as int64) as sb_ref_minus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(1)],",")[OFFSET(0)] as int64) as sb_alt_plus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(1)],",")[OFFSET(1)] as int64) as sb_alt_minus, 
call_AD,
cast(SPLIT(call_AD,",")[OFFSET(0)] as int64) as ref_ad, 
cast(SPLIT(call_AD,",")[OFFSET(1)] as int64) as ad
from position2

union all

select location, sample_id, 
SPLIT(minimize(ref, SPLIT(alt,",")[OFFSET(1)]))[OFFSET(0)] as ref,
SPLIT(minimize(ref, SPLIT(alt,",")[OFFSET(1)]))[OFFSET(1)] as allele,
2 as allele_pos, call_GT,
as_raw_mq,
cast(SPLIT(as_raw_mq,"|")[OFFSET(2)] as int64) raw_mq,
as_raw_mqranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_mqranksum,",")[SAFE_OFFSET(1)] as float64) * 10.0 as int64) as raw_mqranksum_x_10,  # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(1) for the ALT.  We also have some NaNs 
as_qualapprox,
cast(SPLIT(as_qualapprox,"|")[OFFSET(1)] as int64) as qual,  # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(1) for the ALT.  We also have some NaNs 
as_raw_readposranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_readposranksum,",")[SAFE_OFFSET(1)] as float64) * 10.0 as int64) as raw_readposranksum_x_10,  # KC:11-06-20 it seems the leading | has been stripped off, so it's offset(1) for the ALT.  We also have some aNs 
as_sb_table,
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(0)],",")[OFFSET(0)] as int64) as sb_ref_plus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(0)],",")[OFFSET(1)] as int64) as sb_ref_minus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(2)],",")[OFFSET(0)] as int64) as sb_alt_plus, 
cast(SPLIT(SPLIT(as_sb_table,"|")[OFFSET(2)],",")[OFFSET(1)] as int64) as sb_alt_minus, 
call_AD,
cast(SPLIT(call_AD,",")[OFFSET(0)] as int64) as ref_ad, 
cast(SPLIT(call_AD,",")[OFFSET(2)] as int64) as ad
from position2;    