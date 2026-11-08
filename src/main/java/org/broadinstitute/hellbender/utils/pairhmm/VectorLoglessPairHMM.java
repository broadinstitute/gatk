package org.broadinstitute.hellbender.utils.pairhmm;

import com.fulcrumgenomics.fgkl.pairhmm.FgklPairHmm;
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

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for performing the pair HMM for global alignment using the fgkl native library, which selects the fastest
 * SIMD kernel available on the running CPU and computes on the calling thread.
 */
public final class VectorLoglessPairHMM extends LoglessPairHMM {

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
     * @param args              arguments to the native implementation
     * @throws UserException.HardwareFeatureException if no native library is available for this platform
     */
    public VectorLoglessPairHMM(final PairHMMNativeArguments args) throws UserException.HardwareFeatureException {
        final FgklPairHmm fgklPairHmm = new FgklPairHmm();
        if (!fgklPairHmm.load(null)) {
            throw new UserException.HardwareFeatureException("The native PairHMM library is not available on this platform.");
        }
        fgklPairHmm.initialize(args);
        logger.info("Using the native PairHMM with the " + fgklPairHmm.backend() + " backend");
        pairHmm = fgklPairHmm;
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
    public void computeLog10Likelihoods(final LikelihoodMatrix<GATKRead, Haplotype> logLikelihoods,
                                        final List<GATKRead> processedReads, final PairHMMInputScoreImputator inputScoreImputator) {
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
            final PairHMMInputScoreImputation inputScoreImputation = inputScoreImputator.impute(read);
            readDataArray[idx] = new ReadDataHolder();
            readDataArray[idx].readBases = read.getBases();
            readDataArray[idx].readQuals = read.getBaseQualities();
            readDataArray[idx].insertionGOP = inputScoreImputation.insOpenPenalties();
            readDataArray[idx].deletionGOP = inputScoreImputation.delOpenPenalties();
            readDataArray[idx].overallGCP = inputScoreImputation.gapContinuationPenalties();
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
                writeToResultsFileIfApplicable(readDataArray[r].readBases, readDataArray[r].readQuals, readDataArray[r].insertionGOP, readDataArray[r].deletionGOP, readDataArray[r].overallGCP, haplotype.getBases(), mLogLikelihoodArray[readIdx + idxInsideHaplotypeList]);
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
