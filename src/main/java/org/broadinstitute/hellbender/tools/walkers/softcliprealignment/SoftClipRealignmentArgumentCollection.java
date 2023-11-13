package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;

public class SoftClipRealignmentArgumentCollection {
    public static final int DEFAULT_MIN_SEED_LENGTH = 14;
    public static final double DEFAULT_DROP_RATIO = 0.2;
    public static final double DEFAULT_SEED_SPLIT_FACTOR = 0.5;
    public static final int DEFAULT_MAX_REASONABLE_FRAGMENT_LENGTH = 100000;
    public static final double DEFAULT_MIN_ALIGNER_SCORE_DIFFERENCE_PER_BASE = 0.2;
    public static final double DEFAULT_MIN_MISMATCH_DIFFERENCE_PER_BASE = 0.02;
    public static final int DEFAULT_NUM_REGULAR_CONTIGS = 250000;


    @Argument(doc="BWA index image path generated with BwaMemIndexImageCreator.",
            fullName = BwaArgumentCollection.BWA_MEM_INDEX_IMAGE_FULL_NAME)
    public GATKPath indexImage;

    @Argument(doc="Number of reads to hold in the buffer before aligning as a batch. This should be sufficiently " +
            "large that the aligner can accurately estimate the insert size distribution.",
            fullName = "buffer-size",
            optional = true,
            minValue = 100)
    public int bufferSize = 40000;

    @Argument(doc="Number of bwa threads.",
            fullName = "bwa-threads",
            optional = true,
            minValue = 1)
    public int bwaThreads = 2;

}