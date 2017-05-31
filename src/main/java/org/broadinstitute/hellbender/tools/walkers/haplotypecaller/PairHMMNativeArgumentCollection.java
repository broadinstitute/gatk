package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.barclay.argparser.Argument;

/**
 * Arguments for native PairHMM implementations
 */
public class PairHMMNativeArgumentCollection {

    @Argument(fullName = "nativePairHmmThreads", shortName = "threads", doc="How many threads should a native pairHMM implementation use", optional = true)
    private int pairHmmNativeThreads = 1;

    @Argument(fullName = "useDoublePrecision", shortName = "useDoublePrecision", doc="use double precision in the native pairHmm. " +
            "This is slower but matches the java implementation better", optional = true)
    private boolean useDoublePrecision = false;

    public PairHMMNativeArguments getPairHMMArgs(){
        final PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.maxNumberOfThreads = pairHmmNativeThreads;
        args.useDoublePrecision = useDoublePrecision;
        return args;
    }

}
