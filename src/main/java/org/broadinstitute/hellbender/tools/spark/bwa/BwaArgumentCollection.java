package org.broadinstitute.hellbender.tools.spark.bwa;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;

/**
 * Arguments for BWA.
 */
public final class BwaArgumentCollection implements ArgumentCollectionDefinition {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "the number of base pairs to send in a batch to BWA", shortName = "K",
            fullName = "fixedChunkSize", optional = true)
    public int fixedChunkSize = 100000;

    @Argument(doc = "the number of threads", shortName = "t",
            fullName = "threads", optional = true)
    public int numThreads = 1;
}
