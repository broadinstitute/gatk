package org.broadinstitute.hellbender.utils.pairhmm;

import com.intel.gkl.pairhmm.IntelPairHmm;
import com.intel.gkl.pairhmm.IntelPairHmmOMP;
import com.intel.gkl.pairhmm.IntelPairHmmFpga;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeBinding;
import org.broadinstitute.gatk.nativebindings.pairhmm.ReadDataHolder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for performing the pair HMM for local alignment using AVX instructions contained in a native shared library.
 */
public final class VectorLoglessPairHMM extends LoglessPairHMM {

    /**
     * Type for implementation of VectorLoglessPairHMM
     */
    public enum Implementation {
        /**
         * AVX-accelerated version of PairHMM
         */
        AVX,
        /**
         * OpenMP multi-threaded AVX-accelerated version of PairHMM
         */
        OMP,
        /**
         * FPGA-accelerated version of PairHMM
         */
        FPGA
    }

    private static final Logger logger = LogManager.getLogger(VectorLoglessPairHMM.class);
    private long threadLocalSetupTimeDiff = 0;
    private long pairHMMSetupTime = 0;

    private final PairHMMNativeBinding pairHmm;

    //Hold the mapping between haplotype and index in the list of Haplotypes passed to initialize
    //Use this mapping in computeLikelihoods to find the likelihood value corresponding to a given Haplotype
    private final Map<Haplotype, Integer> haplotypeToHaplotypeListIdxMap = new LinkedHashMap<>();
    private HaplotypeDataHolder[] mHaplotypeDataArray;

    /**
     * Create a VectorLoglessPairHMM
     *
     * @param implementation    which implementation to use (AVX or OMP)
     * @param args              arguments to the native GKL implementation
     */
    public VectorLoglessPairHMM(Implementation implementation, PairHMMNativeArguments args) throws UserException.HardwareFeatureException {
        final boolean isSupported;

        switch (implementation) {
            case AVX:
                pairHmm = new IntelPairHmm();
                isSupported = pairHmm.load(null);
                if (!isSupported) {
                    throw new UserException.HardwareFeatureException("Machine does not support AVX PairHMM.");
                }
                break;

            case OMP:
                pairHmm = new IntelPairHmmOMP();
                isSupported = pairHmm.load(null);
                if (!isSupported) {
                    throw new UserException.HardwareFeatureException("Machine does not support OpenMP AVX PairHMM.");
                }
                break;

            case FPGA:
                pairHmm = new IntelPairHmmFpga();
                isSupported = pairHmm.load(null);
                if (!isSupported) {
                    throw new UserException.HardwareFeatureException("Machine does not support FPGA PairHMM.");
                }
                break;

            default:
                throw new UserException.HardwareFeatureException("Unknown PairHMM implementation.");
        }

        pairHmm.initialize(args);
    }


    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize(final List<Haplotype> haplotypes, final Map<String, List<GATKRead>> perSampleReadList,
                           final int readMaxLength, final int haplotypeMaxLength) {
        // do not need to call super.initialize()
        int numHaplotypes = haplotypes.size();
        mHaplotypeDataArray = new HaplotypeDataHolder[numHaplotypes];
        int idx = 0;
        haplotypeToHaplotypeListIdxMap.clear();
        for (final Haplotype currHaplotype : haplotypes) {
            mHaplotypeDataArray[idx] = new HaplotypeDataHolder();
            mHaplotypeDataArray[idx].haplotypeBases = currHaplotype.getBases();
            haplotypeToHaplotypeListIdxMap.put(currHaplotype, idx);
            ++idx;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void computeLog10Likelihoods(final LikelihoodMatrix<Haplotype> logLikelihoods,
                                        final List<GATKRead> processedReads,
                                        final Map<GATKRead, byte[]> gcp) {
        if (processedReads.isEmpty()) {
            return;
        }
        if (doProfiling) {
            startTime = System.nanoTime();
        }
        int readListSize = processedReads.size();
        int numHaplotypes = logLikelihoods.numberOfAlleles();
        ReadDataHolder[] readDataArray = new ReadDataHolder[readListSize];
        int idx = 0;
        for (GATKRead read : processedReads) {
            readDataArray[idx] = new ReadDataHolder();
            readDataArray[idx].readBases = read.getBases();
            readDataArray[idx].readQuals = read.getBaseQualities();
            readDataArray[idx].insertionGOP = ReadUtils.getBaseInsertionQualities(read);
            readDataArray[idx].deletionGOP = ReadUtils.getBaseDeletionQualities(read);
            readDataArray[idx].overallGCP = gcp.get(read);
            ++idx;
        }

        mLogLikelihoodArray = new double[readListSize * numHaplotypes];      //to store results
        if (doProfiling) {
            threadLocalSetupTimeDiff = (System.nanoTime() - startTime);
        }
        //for(reads)
        //   for(haplotypes)
        //       compute_full_prob()
        pairHmm.computeLikelihoods(readDataArray, mHaplotypeDataArray, mLogLikelihoodArray);

        int readIdx = 0;
        for (int r = 0; r < readListSize; r++) {
            int hapIdx = 0;
            for (final Haplotype haplotype : logLikelihoods.alleles()) {

                //Since the order of haplotypes in the List<Haplotype> and alleleHaplotypeMap is different,
                //get idx of current haplotype in the list and use this idx to get the right likelihoodValue
                final int idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.get(haplotype);
                logLikelihoods.set(hapIdx, r, mLogLikelihoodArray[readIdx + idxInsideHaplotypeList]);
                ++hapIdx;
            }
            readIdx += numHaplotypes;
        }
        if (doProfiling) {
            threadLocalPairHMMComputeTimeDiff = (System.nanoTime() - startTime);
            pairHMMComputeTime += threadLocalPairHMMComputeTimeDiff;
            pairHMMSetupTime += threadLocalSetupTimeDiff;
        }
    }


    @Override
    public void close() {
        pairHmm.done();
        if (doProfiling)
            logger.info("Time spent in setup for JNI call : " + (pairHMMSetupTime * 1e-9));
        super.close();
    }
}
