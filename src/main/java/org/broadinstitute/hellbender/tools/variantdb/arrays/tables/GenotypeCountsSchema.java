package org.broadinstitute.hellbender.tools.variantdb.arrays.tables;

import java.util.Arrays;
import java.util.List;

public class GenotypeCountsSchema {

    public static final String PROBE_ID = "probe_id";
    public static final String HOM_REF_COUNT = "hom_ref";
    public static final String HET_COUNT = "het";
    public static final String HOM_VAR_COUNT = "hom_var";
    public static final String HET_1_2_COUNT = "het_1_2";
    public static final String HOM_VAR_2_2_COUNT = "hom_var_2_2";
    public static final String NO_CALL_COUNT = "no_call";


    public static final List<String> GENOTYPE_COUNTS_FIELDS = Arrays.asList(PROBE_ID, HOM_REF_COUNT, HET_COUNT, HOM_VAR_COUNT, HET_1_2_COUNT, HOM_VAR_2_2_COUNT, NO_CALL_COUNT);
}
