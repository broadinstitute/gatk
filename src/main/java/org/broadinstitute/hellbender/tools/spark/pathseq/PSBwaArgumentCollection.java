package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public final class PSBwaArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "URI of the pathogen BWA image file for alignment. This must be distributed to local disk on each node.",
            fullName = "pathogenBwaImage")
    public String bwaImage;

    @Argument(doc = "Reference URI corresponding to the pathogen image file.",
            fullName = "pathogenFasta")
    public String referencePath;

    @Argument(doc = "Minimum seed length for the pathogen BWA alignment.",
            fullName = "pathogenMinSeedLength",
            minValue = 1,
            minRecommendedValue = 11,
            optional = true)
    public int seedLength = 19;

    @Argument(doc = "Number of BWA threads for pathogen alignment. In Spark, this is the number of threads used per task.",
            fullName = "bwaThreads",
            minValue = 1,
            optional = true)
    public int bwaThreads = 1;

    @Argument(doc = "Maximum number of alternate pathogen hits.",
            fullName = "maxAlternateHits",
            minValue = 0,
            optional = true)
    public int maxAlternateHits = 5000;

    @Argument(doc = "Minimum score threshold for pathogen alignment output.",
            fullName = "hostScoreThreshold",
            minValue = 0,
            optional = true)
    public int scoreThreshold = 30;

}