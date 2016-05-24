package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Constains utilities for tests that need to create read-likelihoods.
 */
public final class ReadLikelihoodsUnitTester {

    public static ReadLikelihoods<Allele> readLikelihoods(final int alleleCount, final int[] readCount) {
        final int sampleCount = readCount.length;
        final AlleleList<Allele> alleleList = AlleleListUnitTester.alleleList(alleleCount, 100, true);
        final SampleList sampleList = SampleListUnitTester.sampleList(sampleCount);
        final Map<String,List<GATKRead>> sampleToReads = new LinkedHashMap<>(sampleCount);
        for (int i = 0; i < sampleCount; i++) {
            sampleToReads.put(sampleList.getSample(i),readList(i,readCount[i]));
        }
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList,alleleList, sampleToReads);
        for (int s = 0; s < sampleCount; s++) {
            final LikelihoodMatrix<Allele> sampleLikelihoods = likelihoods.sampleMatrix(s);
            for (int a = 0; a < alleleCount; a++)
                for (int r = 0; r < readCount[s]; r++)
                    sampleLikelihoods.set(a, r, testLikelihood(s, a, r));
        }
        return likelihoods;
    }

    /**
     * produces a test likelihood depending on the sample, read and allele index.
     */
    private static double testLikelihood(final int sampleIndex, final int alleleIndex, final int readIndex) {
        return - Math.abs(3 * (sampleIndex + 1) + 7 * (alleleIndex + 1) + 11 * (readIndex + 1));
    }


    private static SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(10, 0, 1000);


    static List<GATKRead> readList(final int sampleIndex, final int readCount) {
        final List<GATKRead> reads = new ArrayList<>(readCount);
        int readIndex = 0;
        for (int j = 0; j < readCount; j++)
            reads.add(ArtificialReadUtils.createArtificialRead(SAM_HEADER, "READ_" + sampleIndex + "_" + (readIndex++), 1, 1, 100));
        return reads;
    }

    /**
     * Creates a sampleToReads map given the sample list and the required read counts.
     * @param sampleList the target sample-list.
     * @param readCounts the target read-counts.
     * @return never {@code null}.
     */
    public static Map<String,List<GATKRead>> sampleToReads(final SampleList sampleList, final int[] readCounts) {
        final Map<String,List<GATKRead>> result = new LinkedHashMap<>(sampleList.numberOfSamples());
        int readIndex = 0;
        for (int i = 0; i < sampleList.numberOfSamples(); i++) {
            final int readCount = readCounts[i];
            final String sample = sampleList.getSample(i);
            final List<GATKRead> records = new ArrayList<>(readCount);
            for (int j = 0; j < readCount; j++)
                records.add(ArtificialReadUtils.createArtificialRead(SAM_HEADER,"READ_" + (readIndex++),1,1,100));
            result.put(sample,records);
        }
        return result;
    }

}