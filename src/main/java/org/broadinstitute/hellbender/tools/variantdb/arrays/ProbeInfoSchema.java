package org.broadinstitute.hellbender.tools.variantdb.arrays;

import java.sql.Ref;
import java.util.Arrays;
import java.util.List;

public class ProbeInfoSchema {

    public static final String PROBE_ID = "ProbeId";
    public static final String NAME = "Name";
    public static final String REF_BUILD = "GenomeBuild";
    public static final String CONTIG = "Chr";
    public static final String POSITION = "Position";
    public static final String REF = "Ref";
    public static final String ALLELE_A = "AlleleA";
    public static final String ALLELE_B = "AlleleB";
    public static final String FLAG = "build37Flag";

    public static final List<String> PROBE_INFO_FIELDS = Arrays.asList(PROBE_ID, NAME, CONTIG, POSITION, REF, ALLELE_A, ALLELE_B, FLAG);
    public static final List<String> PROBE_INFO_FOR_INGEST_FIELDS = Arrays.asList(PROBE_ID, NAME, REF_BUILD, REF, ALLELE_A, ALLELE_B);

}
