package org.broadinstitute.hellbender.tools.funcotator;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class FuncotatorConstants {

    /**
     * A character representing any base as used in VCF files and IUPAC base standards.
     * NOTE: When HTSJDK with the fix for issue #1226 (https://github.com/samtools/htsjdk/issues/1226) is released
     *       we will want to remove this constant.
     */
    public static final char MASKED_ANY_BASE = 'N';
    /**
     * A {@link String} representing any base as used in VCF files and IUPAC base standards.
     * NOTE: When HTSJDK with the fix for issue #1226 (https://github.com/samtools/htsjdk/issues/1226) is released
     *       we will want to remove this constant.
     */
    public static final String MASKED_ANY_BASE_STRING = String.valueOf(MASKED_ANY_BASE);

    /**
     * Value to insert into fields when the annotation is unspecified, but still required.
     */
    public static final String UNKNOWN_VALUE_STRING = "__UNKNOWN__";

    /**
     *
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

    /** Name of the header line key which indicates that a file has already been annotated. */
    public static final String FUNCOTATOR_VERSION_VCF_HEADERLINE_KEY = "Funcotator Version";


}
