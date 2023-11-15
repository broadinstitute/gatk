package org.broadinstitute.hellbender.tools.walkers.softcliprealignment;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;

public class SubsettingRealignmentArgumentCollection {

    public static final String BWA_IMAGE_LONG_NAME = BwaArgumentCollection.BWA_MEM_INDEX_IMAGE_FULL_NAME;
    public static final String BUFFER_SIZE_LONG_NAME = "buffer-size";
    public static final String BWA_THREADS_LONG_NAME = "bwa-threads";

    @Argument(doc="BWA index image path generated with BwaMemIndexImageCreator.",
            fullName = BwaArgumentCollection.BWA_MEM_INDEX_IMAGE_FULL_NAME)
    public GATKPath indexImage;

    @Argument(doc="Number of reads to hold in the buffer before aligning as a batch. This should be sufficiently " +
            "large that the aligner can accurately estimate the insert size distribution.",
            fullName = BUFFER_SIZE_LONG_NAME,
            optional = true,
            minValue = 100)
    public int bufferSize = 40000;

    @Argument(doc="Number of bwa threads.",
            fullName = BWA_THREADS_LONG_NAME,
            optional = true,
            minValue = 1)
    public int bwaThreads = 2;

}