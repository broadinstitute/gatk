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
