package com.intel.gkl.pairhmm;

/**
 * Provides a native PairHMM implementation accelerated using Intel FPGAs
 */
public final class IntelPairHmmFpga extends IntelPairHmm {
    private static final String NATIVE_LIBRARY_NAME = "gkl_pairhmm_fpga";

    public IntelPairHmmFpga() {
        setNativeLibraryName(NATIVE_LIBRARY_NAME);
        useFpga = true;
    }
}
