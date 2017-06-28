package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public final class PSBwaArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "URI of the BWA index image file",
            fullName = "bwamemIndexImage")
    public String bwaImage;

    @Argument(doc = "Reference URI corresponding to the image file",
            fullName = "referencePath")
    public String referencePath;

    @Argument(doc = "Minimum seed length for the BWA alignment.",
            fullName = "minSeedLength",
            minValue = 1,
            minRecommendedValue = 11,
            optional = true)
    public int seedLength = 19;

    @Argument(doc = "Number of BWA threads. In Spark, this is the number of threads used per task.",
            fullName = "bwaThreads",
            minValue = 1,
            optional = true)
    public int bwaThreads = 1;

    @Argument(doc = "Maximum number of alternate hits",
            fullName = "maxAlternateHits",
            minValue = 0,
            optional = true)
    public int maxAlternateHits = 5000;

    @Argument(doc = "Minimum score threshold for alignment output",
            fullName = "scoreThreshold",
            minValue = 0,
            optional = true)
    public int scoreThreshold = 30;

}