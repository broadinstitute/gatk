package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pdhmm.PDHMMNativeArguments;

/**
 * Arguments for native PairHMM implementations
 */
public class PDPairHMMNativeArgumentCollection {

    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-hmm-threads", doc="How many threads should a native DRAGEN pd-pairHMM implementation use", optional = true)
    private int pairHmmNativeThreads = 8;

    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-avx-level", doc="AVX level instruction for gkl pd-pair-hmm impelmentation to use [SCALAR, AVX2, AVX512].", optional = true)
    private PDHMMNativeArguments.AVXLevel avxLevel = PDHMMNativeArguments.AVXLevel.FASTEST_AVAILABLE;

    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-memory-limit", doc="Maximum limit for memory (in MB) to the GKL for internal uses. (NOTE: this memory is not manage by Java and is not subject to -Xmx limits)", optional = true)
    private int memoryLimit = 1024;

    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-open-mp-settings", doc="Whether to enable use of Open MP library within the GKL.", optional = true)
    private PDHMMNativeArguments.OpenMPSetting openMP = PDHMMNativeArguments.OpenMPSetting.FASTEST_AVAILABLE;

    public PDHMMNativeArguments getPDPairHMMArgs(){
        final PDHMMNativeArguments args = new PDHMMNativeArguments();
        args.maxNumberOfThreads = pairHmmNativeThreads;
        args.avxLevel = avxLevel;
        args.setMaxMemoryInMB(memoryLimit);
        args.openMPSetting = openMP;
        return args;
    }

}
