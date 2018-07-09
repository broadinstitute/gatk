package org.broadinstitute.hellbender.utils.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Dummy Read-Likelihood container that computes partial likelihoods based on the read pileups for snps.
 *
 * Note: this class uses FastUtil collections for speed.
 */
public class UnfilledReadsLikelihoods<A extends Allele> extends ReadLikelihoods<A> {

    private Map<String, List<PileupElement>> stratifiedPileups;

    /**
     * Constructs a new read-likelihood collection.
     * <p>
     * <p>
     * The initial likelihoods for all allele-read combinations are
     * 0.
     * </p>
     *
     * @param samples all supported samples in the collection.
     * @param alleles all supported alleles in the collection.
     * @param reads   reads stratified per sample.
     * @throws IllegalArgumentException if any of {@code allele}, {@code samples}
     *                                  or {@code reads} is {@code null},
     *                                  or if they contain null values.
     */
    public UnfilledReadsLikelihoods(SampleList samples, AlleleList<A> alleles, Map<String, List<GATKRead>> reads) {
        super(samples, alleles, reads);
    }

    // Internally used constructor.
    @SuppressWarnings({"unchecked", "rawtypes"})
    private UnfilledReadsLikelihoods(final AlleleList alleles,
                            final SampleList samples,
                            final GATKRead[][] readsBySampleIndex,
                            final Object2IntMap<GATKRead>[] readIndex,
                            final double[][][] values) {
       super(alleles, samples, readsBySampleIndex, readIndex, values);
    }

    /**
     * Is this container expected to have the per-allele likelihoods calculations filled in.
     * This will always be false for this class because we can't compute the likelihoods without re-genotyping the reads
     */
    @Override
    public boolean hasFilledLikelihoods() {
        return false;
    }

    /**
     * Collect a map stratified per-sample of the base pileups at the provided Location
     *
     * @param loc reference location to construct pileups for
     * @return
     */
    public Map<String, List<PileupElement>> getStratifiedPileups(final Locatable loc) {
        if (stratifiedPileups != null) {
            return stratifiedPileups;
        }
        Map<String, List<PileupElement>> pileupMap = new HashMap<>();
        for (int i = 0; i < samples.numberOfSamples(); i++) {
            pileupMap.put(samples.getSample(i), ReadPileup.locToReadsPileup(sampleReads(i), loc));
        }
        stratifiedPileups = Collections.unmodifiableMap(pileupMap);
        return stratifiedPileups;
    }

    /**
     * Create an independent copy of this read-likelihoods collection
     */
    @Override
    ReadLikelihoods<A> copy() {

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        final double[][][] newLikelihoodValues = new double[sampleCount][alleleCount][];

        @SuppressWarnings({"unchecked", "rawtypes"})
        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];

        for (int s = 0; s < sampleCount; s++) {
            newReadsBySampleIndex[s] = readsBySampleIndex[s].clone();
            for (int a = 0; a < alleleCount; a++) {
                newLikelihoodValues[s][a] = valuesBySampleIndex[s][a].clone();
            }
        }

        // Finally we create the new read-likelihood
        return new UnfilledReadsLikelihoods<A>(
                alleles,
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues);
    }

    // Methods Which Modify Reads that must be turned off
    @Override
    public void changeReads(final Map<GATKRead, GATKRead> readRealignments) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }

    @Override
    public void filterPoorlyModeledReads(final double maximumErrorPerBase) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }

    @Override
    public void addReads(final Map<String,List<GATKRead>> readsBySample, final double initialLikelihood) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }

    @Override
    public void contaminationDownsampling(final Map<String, Double> perSampleDownsamplingFraction) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }

    @Override
    public void filterToOnlyOverlappingReads(final SimpleInterval location) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }
}
