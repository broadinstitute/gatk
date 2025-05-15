package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.gatk.nativebindings.pdhmm.PDHMMNativeArguments;

/**
 * Arguments for native PairHMM implementations
 */
public class PDPairHMMNativeArgumentCollection {

    public static final int THREADS = 8;
    public static final PDHMMNativeArguments.AVXLevel AVX_LEVEL = PDHMMNativeArguments.AVXLevel.FASTEST_AVAILABLE;
    public static final int MEMLIMIT = 1024;
    public static final PDHMMNativeArguments.OpenMPSetting OPEN_MP_SETTING = PDHMMNativeArguments.OpenMPSetting.FASTEST_AVAILABLE;
    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-hmm-threads", doc="How many threads should a native DRAGEN pd-pairHMM implementation use", optional = true)
    private int pairHmmNativeThreads = THREADS;

    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-avx-level", doc="AVX level instruction for gkl pd-pair-hmm impelmentation to use [SCALAR, AVX2, AVX512].", optional = true)
    private PDHMMNativeArguments.AVXLevel avxLevel = AVX_LEVEL;

    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-memory-limit", doc="Maximum limit for memory (in MB) to the GKL for internal uses. (NOTE: this memory is not manage by Java and is not subject to -Xmx limits)", optional = true)
    private int memoryLimit = MEMLIMIT;

    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-open-mp-settings", doc="Whether to enable use of Open MP library within the GKL.", optional = true)
    private PDHMMNativeArguments.OpenMPSetting openMP = OPEN_MP_SETTING;

    // NOTE: the PDHMMNativeArguments class is in the GKL code and thus we must work around it here.
    public PDHMMNativeArguments getPDPairHMMArgs(){
        final PDHMMNativeArguments args = new PDHMMNativeArguments();
        args.maxNumberOfThreads = pairHmmNativeThreads;
        args.avxLevel = avxLevel;
        args.setMaxMemoryInMB(memoryLimit);
        args.openMPSetting = openMP;
        return args;
    }
    // The args below are defaults provided for the GKL. The GKL arguments class does not have arguments
    public static PDHMMNativeArguments getDefaultPDPairHMMArgs() {
        final PDHMMNativeArguments args = new PDHMMNativeArguments();
        args.maxNumberOfThreads = THREADS;
        args.avxLevel = AVX_LEVEL;
        args.setMaxMemoryInMB(MEMLIMIT);
        args.openMPSetting = OPEN_MP_SETTING;
        return args;
    }
}
