package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.DeprecatedFeature;
import org.broadinstitute.gatk.nativebindings.pdhmm.PDHMMNativeArguments;

/**
 * Arguments for native PDPairHMM implementations. The native implementation runs on the calling thread and selects
 * its own kernel and memory, so every argument here is accepted for compatibility and has no effect.
 */
public class PDPairHMMNativeArgumentCollection {
    private static final Logger logger = LogManager.getLogger(PDPairHMMNativeArgumentCollection.class);

    public static final int THREADS = 1;
    public static final PDHMMNativeArguments.AVXLevel AVX_LEVEL = PDHMMNativeArguments.AVXLevel.FASTEST_AVAILABLE;
    public static final int MEMLIMIT = 1024;
    public static final PDHMMNativeArguments.OpenMPSetting OPEN_MP_SETTING = PDHMMNativeArguments.OpenMPSetting.FASTEST_AVAILABLE;

    @DeprecatedFeature(detail = "The native PDPairHMM runs on the calling thread, so this argument has no effect")
    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-hmm-threads", doc="Ignored: the native DRAGEN pd-pairHMM runs on the calling thread", optional = true)
    private int pairHmmNativeThreads = THREADS;

    @DeprecatedFeature(detail = "The native PDPairHMM selects the fastest kernel for the running CPU, so this argument has no effect")
    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-avx-level", doc="Ignored: the native DRAGEN pd-pairHMM selects the fastest kernel for the running CPU", optional = true)
    private PDHMMNativeArguments.AVXLevel avxLevel = AVX_LEVEL;

    @DeprecatedFeature(detail = "The native PDPairHMM manages its own memory, so this argument has no effect")
    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-memory-limit", doc="Ignored: the native DRAGEN pd-pairHMM manages its own memory", optional = true)
    private int memoryLimit = MEMLIMIT;

    @DeprecatedFeature(detail = "The native PDPairHMM runs on the calling thread, so this argument has no effect")
    @Advanced
    @Argument(fullName = "native-dragen-pd-pair-open-mp-settings", doc="Ignored: the native DRAGEN pd-pairHMM runs on the calling thread", optional = true)
    private PDHMMNativeArguments.OpenMPSetting openMP = OPEN_MP_SETTING;

    public PDHMMNativeArguments getPDPairHMMArgs(){
        if (pairHmmNativeThreads != THREADS || avxLevel != AVX_LEVEL || memoryLimit != MEMLIMIT || openMP != OPEN_MP_SETTING) {
            logger.warn("The --native-dragen-pd-pair-* arguments are deprecated and have no effect: the native PDPairHMM runs on the calling thread and selects its own kernel");
        }
        final PDHMMNativeArguments args = new PDHMMNativeArguments();
        args.maxNumberOfThreads = pairHmmNativeThreads;
        args.avxLevel = avxLevel;
        args.setMaxMemoryInMB(memoryLimit);
        args.openMPSetting = openMP;
        return args;
    }

    // The gatk-native-bindings arguments class has no defaults of its own
    public static PDHMMNativeArguments getDefaultPDPairHMMArgs() {
        final PDHMMNativeArguments args = new PDHMMNativeArguments();
        args.maxNumberOfThreads = THREADS;
        args.avxLevel = AVX_LEVEL;
        args.setMaxMemoryInMB(MEMLIMIT);
        args.openMPSetting = OPEN_MP_SETTING;
        return args;
    }
}
