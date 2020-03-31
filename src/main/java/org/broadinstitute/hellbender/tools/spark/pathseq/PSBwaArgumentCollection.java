package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public final class PSBwaArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String MICROBE_BWA_IMAGE_LONG_NAME = "microbe-bwa-image";
    public static final String MICROBE_REF_DICT_LONG_NAME = "microbe-dict";
    public static final String MICROBE_MIN_SEED_LENGTH_LONG_NAME = "microbe-min-seed-length";
    public static final String MAX_ALT_HITS_LONG_NAME = "max-alternate-hits";
    public static final String SCORE_THRESHOLD_LONG_NAME = "bwa-score-threshold";

    @Argument(doc = "Microbe reference BWA index image file generated using BwaMemIndexImageCreator. If running on a Spark cluster, this must be distributed to local disk on each node.",
            fullName = MICROBE_BWA_IMAGE_LONG_NAME)
    public String bwaImage;

    @Argument(fullName = MICROBE_REF_DICT_LONG_NAME,
            doc = "Use the given sequence dictionary as the microbe sequence dictionary.  Must be a .dict file.")
    public String microbeDictionary = null;

    /**
     * This parameter controls the sensitivity of the BWA-MEM aligner. Smaller values result in more alignments at the
     * expense of computation time.
     */
    @Argument(doc = "Minimum BWA-MEM seed length for the microbe alignment",
            fullName = MICROBE_MIN_SEED_LENGTH_LONG_NAME,
            minValue = 1,
            minRecommendedValue = 11,
            optional = true)
    public int seedLength = 19;

    /**
     * The maximum number of alternate alignments for each read, i.e. the alignments appearing in the XA tag.
     */
    @Argument(doc = "Maximum number of alternate microbe alignments",
            fullName = MAX_ALT_HITS_LONG_NAME,
            minValue = 0,
            optional = true)
    public int maxAlternateHits = 5000;

    /**
     * This parameter controls the minimum quality of the BWA alignments for the output.
     */
    @Argument(doc = "Minimum score threshold for microbe alignments",
            fullName = SCORE_THRESHOLD_LONG_NAME,
            minValue = 0,
            optional = true)
    public int scoreThreshold = 30;

    public final int bwaThreads = 1;

}