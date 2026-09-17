package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.DeprecatedFeature;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;

/**
 * Arguments for native PairHMM implementations
 */
public class PairHMMNativeArgumentCollection {
    private static final Logger logger = LogManager.getLogger(PairHMMNativeArgumentCollection.class);

    @DeprecatedFeature(detail = "The native PairHMM runs on the calling thread, so this argument has no effect")
    @Argument(fullName = "native-pair-hmm-threads", doc="Ignored: the native PairHMM runs on the calling thread", optional = true)
    private int pairHmmNativeThreads = 1;

    @Argument(fullName = "native-pair-hmm-use-double-precision", doc="use double precision in the native pairHmm. " +
            "This is slower but matches the java implementation better", optional = true)
    private boolean useDoublePrecision = false;

    public PairHMMNativeArguments getPairHMMArgs(){
        if (pairHmmNativeThreads != 1) {
            logger.warn("--native-pair-hmm-threads is deprecated and has no effect: the native PairHMM runs on the calling thread");
        }
        final PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.maxNumberOfThreads = pairHmmNativeThreads;
        args.useDoublePrecision = useDoublePrecision;
        return args;
    }

}
