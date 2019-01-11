package org.broadinstitute.hellbender.tools.spark.bwa;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.BwaMemIndexImageCreator;

import java.io.Serializable;

/**
 * A collection of the arguments that are used for BWA.
 */
public class BwaArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String SINGLE_END_ALIGNMENT_FULL_NAME = "single-end-alignment";
    public static final String SINGLE_END_ALIGNMENT_SHORT_NAME = "se";
    public static final String BWA_MEM_INDEX_IMAGE_FULL_NAME = "bwa-mem-index-image";
    public static final String BWA_MEM_INDEX_IMAGE_SHORT_NAME = "image";

    /**
     * The BWA-MEM index image file name that you've distributed to each executor. The image file can be generated using
     * {@link BwaMemIndexImageCreator}. If this argument is not specified, the default behavior is to look for a
     * file whose name is the FASTA reference file with a <i>.img</i> suffix.
     */
    @Argument(doc = "The BWA-MEM index image file name that you've distributed to each executor",
            fullName = BWA_MEM_INDEX_IMAGE_FULL_NAME,
            shortName = BWA_MEM_INDEX_IMAGE_SHORT_NAME,
            optional = true)
    public String indexImageFile;

    /**
     * Run single-end instead of paired-end alignment.
     */
    @Argument(doc = "Run single-end instead of paired-end alignment",
            fullName = SINGLE_END_ALIGNMENT_FULL_NAME,
            shortName = SINGLE_END_ALIGNMENT_SHORT_NAME,
            optional = true)
    public boolean singleEndAlignment = false;
}
