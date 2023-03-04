package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

public class FlowBasedAlignmentArgumentCollection extends FlowBasedArgumentCollection {

    private static final long serialVersionUID = 0;

    public static final String FLOW_LIKELIHOOD_PARALLEL_THREADS_LONG_NAME = "flow-likelihood-parallel-threads";
    public static final String FLOW_LIKELIHOOD_OPTIMIZED_COMP_LONG_NAME = "flow-likelihood-optimized-comp";

    public static final String TRIM_TO_HAPLOTYPE_LONG_NAME = "trim-to-haplotype";
    public static final String EXACT_MATCHING_LONG_NAME = "exact-matching";
    @Advanced
    @Hidden
    @Argument(fullName = FLOW_LIKELIHOOD_PARALLEL_THREADS_LONG_NAME, doc = "Number of threads to parallelize likelihood computation inner (read) loop with", optional=true)
    public int flowLikelihoodParallelThreads = 0;

    @Advanced
    @Hidden
    @Argument(fullName = FLOW_LIKELIHOOD_OPTIMIZED_COMP_LONG_NAME, doc = "Use optimized likelihood computation version. The code is otimized in that it performs fewer log10 calls - which are expensive - by using precomputed values " +
            "for common probability values", optional=true)
    public boolean flowLikelihoodOptimizedComp = false;

    @Advanced
    @Hidden
    @Argument(fullName = TRIM_TO_HAPLOTYPE_LONG_NAME, doc = "Should the read be trimmed to the haplotype", optional=true)
    public boolean trimToHaplotype = true;

    @Advanced
    @Hidden
    @Argument(fullName = EXACT_MATCHING_LONG_NAME, doc = "Should exact match of haplotype to read be assumed", optional=true)
    public boolean exactMatching = false;

}
