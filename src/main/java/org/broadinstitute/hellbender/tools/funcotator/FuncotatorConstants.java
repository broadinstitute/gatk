package org.broadinstitute.hellbender.tools.funcotator;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class FuncotatorConstants {
    /**
     * Datasource name to use for Funcotations created from input variants from a VCF.
     */
    public static final String DATASOURCE_NAME_FOR_INPUT_VCFS = "INPUT_VCF";

    /** Name of the Mitochondrial contig for the B37 human reference. */
    public static final String B37_MITOCHONDRIAL_CONTIG_NAME = "MT";

    /** Name of the Mitochondrial contig for the HG19 human reference. */
    public static final String HG19_MITOCHONDRIAL_CONTIG_NAME = "chrM";

    /** Name of the Mitochondrial contig for the HG38 human reference. */
    public static final String HG38_MITOCHONDRIAL_CONTIG_NAME = "chrM";

    /** {@link Set<String>} containing known valid names for Mitochondrial chromosomes. */
    public static final Set<String> MITOCHONDRIAL_CONTIG_NAMES = new HashSet<>(
            Arrays.asList(B37_MITOCHONDRIAL_CONTIG_NAME, HG19_MITOCHONDRIAL_CONTIG_NAME, HG38_MITOCHONDRIAL_CONTIG_NAME)
    );

}
