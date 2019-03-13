package com.intel.gkl.pairhmm;

/**
 * Provides a native PairHMM implementation accelerated for the Intel Architecture using multithreaded OpenMP.
 */
public final class IntelPairHmmOMP extends IntelPairHmm {
    private static final String NATIVE_LIBRARY_NAME = "gkl_pairhmm_omp";

    public IntelPairHmmOMP() {
        setNativeLibraryName(NATIVE_LIBRARY_NAME);
        useOmp = true;
    }
}
