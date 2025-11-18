package org.broadinstitute.hellbender.utils.pairhmm;

import com.intel.gkl.pdhmm.IntelPDHMM;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.SAMUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pdhmm.PDHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pdhmm.PDHMMNativeArguments.AVXLevel;
import org.broadinstitute.gatk.nativebindings.pdhmm.PDHMMNativeArguments.OpenMPSetting;
import org.broadinstitute.gatk.nativebindings.pdhmm.ReadDataHolder;
import org.broadinstitute.gatk.nativebindings.pdhmm.HaplotypeDataHolder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PDPairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for performing the pair HMM for global alignment using AVX instructions
 * contained in a native shared library.
 */
public final class VectorLoglessPairPDHMM extends LoglessPDPairHMM {

    private static final Logger logger = LogManager.getLogger(VectorLoglessPairPDHMM.class);
    private static long pairHMMTime = 0;
    private static long pairHMMSetupTime = 0;
    private static long postProcessingTime = 0;
    private static long totalComps = 0;

    private final IntelPDHMM pairPDHmm;

    /**
     * Create a VectorLoglessPairHMM
     *
     * @param args           arguments to the native GKL implementation
     */
    public VectorLoglessPairPDHMM(final PDHMMNativeArguments args)
            throws UserException.HardwareFeatureException {
        final boolean isSupported;

        // Check if the native library loads (which internally checks for AVX support)
        pairPDHmm = new IntelPDHMM();
        isSupported = pairPDHmm.load(null); // NOTE: the temp dir defaults to Java.IO.File's temp dir
        if (!isSupported) {
            throw new UserException.HardwareFeatureException("Machine does not support OpenMP AVX PairHMM.");
        }

        pairPDHmm.initialize(args);
    }

    public void initialize(final int readMaxLength, final int haplotypeMaxLength) {
        super.initialize(readMaxLength, haplotypeMaxLength);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void computeLog10Likelihoods(final LikelihoodMatrix<GATKRead, PartiallyDeterminedHaplotype> logLikelihoods,
                                        final List<GATKRead> processedReads,
                                        final PairHMMInputScoreImputator inputScoreImputator,
                                        final int rangeForReadOverlapToDeterminedBases) {

        if (processedReads.isEmpty()) {
            return;
        }
        if (doProfiling) {
            startTime = System.nanoTime();
        }
        // (re)initialize the pairHMM only if necessary
        final List<PartiallyDeterminedHaplotype> alleles = logLikelihoods.alleles();
        final int readCount = processedReads.size();
        final int alleleCount = alleles.size();
        totalComps+= (long) readCount * alleleCount;

        ReadDataHolder[] readDataHolders = new ReadDataHolder[readCount];
        for (int readIndex = 0; readIndex < readCount; readIndex++) {
            GATKRead read = processedReads.get(readIndex);
            final PairHMMInputScoreImputation inputScoreImputation = inputScoreImputator.impute(read);
            readDataHolders[readIndex] = new ReadDataHolder();
            readDataHolders[readIndex].readBases = read.getBases();
            readDataHolders[readIndex].readQuals = read.getBaseQualities();
            readDataHolders[readIndex].insertionGOP = inputScoreImputation.insOpenPenalties();
            readDataHolders[readIndex].deletionGOP = inputScoreImputation.delOpenPenalties();
            readDataHolders[readIndex].overallGCP = inputScoreImputation.gapContinuationPenalties();
        }

        HaplotypeDataHolder[] haplotypeDataHolders = new HaplotypeDataHolder[alleleCount];
        for (int alleleIndex = 0; alleleIndex < alleleCount; alleleIndex++) {
            PartiallyDeterminedHaplotype allele = alleles.get(alleleIndex);
            haplotypeDataHolders[alleleIndex] = new HaplotypeDataHolder();
            haplotypeDataHolders[alleleIndex].haplotypeBases = allele.getBases();
            haplotypeDataHolders[alleleIndex].haplotypePDBases = allele.getAlternateBases();
        }

        mLogLikelihoodArray = new double[readCount * alleleCount];

        if (doProfiling) {
            pairHMMSetupTime += (System.nanoTime() - startTime);
            startTime = System.nanoTime();
        }

        pairPDHmm.computeLikelihoods(readDataHolders, haplotypeDataHolders, mLogLikelihoodArray);

        if (doProfiling) {
            pairHMMTime += (System.nanoTime() - startTime);
            startTime = System.nanoTime();
        }

        int readIndex = 0;
        for (final GATKRead read : processedReads) {
            final Locatable unclippedSpan = (Locatable) read.getTransientAttribute(
                    PDPairHMMLikelihoodCalculationEngine.UNCLIPPED_ORIGINAL_SPAN_ATTR);
            for (int a = 0; a < alleleCount; a++) {
                final PartiallyDeterminedHaplotype allele = alleles.get(a);
                final int likelihoodIndex = readIndex * alleleCount + a;
                if (rangeForReadOverlapToDeterminedBases < 0 ||
                    allele.getMaximumExtentOfSiteDeterminedAlleles().overlapsWithMargin(
                        unclippedSpan, rangeForReadOverlapToDeterminedBases + 1)) {
                    // Do nothing, keep the computed likelihood
                } else {
                    mLogLikelihoodArray[likelihoodIndex] = Double.NEGATIVE_INFINITY;
                }
                logLikelihoods.set(a, readIndex, mLogLikelihoodArray[likelihoodIndex]);
                writeToResultsFileIfApplicable(inputScoreImputator, read, allele, mLogLikelihoodArray[likelihoodIndex]);
            }
            readIndex++;
        }

        if(doProfiling) {
            postProcessingTime += (System.nanoTime() - startTime);
        }


    }

    @Override
    public void close() {
        pairPDHmm.done();
        if (doProfiling) {
            logger.info("Time spent in setup for JNI call : " + (pairHMMSetupTime * 1e-9));
            logger.info("Time spent in JNI call : " + (pairHMMTime * 1e-9));
            logger.info("Time spent in post-processing : " + (postProcessingTime * 1e-9));
            logger.info("Total comparisons : " + totalComps);
        }
        super.close();
    }

    private void writeToResultsFileIfApplicable(final PairHMMInputScoreImputator inputScoreImputator, GATKRead read, PartiallyDeterminedHaplotype allele, double lk)
    {
        if(debugOutputStream == null) {
            return;
        }

        if (inputScoreImputator == null || read == null || allele == null) {
            logger.warn("Null input provided to writeToResultsFileIfApplicable");
            return;
        }

        final PairHMMInputScoreImputation inputScoreImputation = inputScoreImputator.impute(read);
        final byte[] readBases = read.getBases();
        final byte[] readQuals = read.getBaseQualities();
        final byte[] readInsQuals = inputScoreImputation.insOpenPenalties();
        final byte[] readDelQuals = inputScoreImputation.delOpenPenalties();
        final byte[] overallGCP = inputScoreImputation.gapContinuationPenalties();

        final byte[] alleleBases = allele.getBases();
        final byte[] allelePDBases = allele.getAlternateBases();

        try {
            debugOutputStream.write("\n" + new String(alleleBases) + "\t" + Arrays.toString(allelePDBases) + "\t" + new String(readBases) + "\t" + SAMUtils.phredToFastq(readQuals) + "\t" + SAMUtils.phredToFastq(readInsQuals) + "\t" + SAMUtils.phredToFastq(readDelQuals) + "\t" + SAMUtils.phredToFastq(overallGCP) + "\t" + String.format("%e",lk));
        } catch (IOException e) {
            throw new UserException("Error writing to specified HMM results output stream", e);
        }

    }
}
