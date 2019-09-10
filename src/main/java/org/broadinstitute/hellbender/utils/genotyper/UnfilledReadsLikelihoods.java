package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.ToDoubleFunction;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotator;

/**
 * Dummy implementation of likelihoods class, used for GATK annotation engine, that doesn't actually contain likelihoods.
 *
 * Several GATK annotations use the sample read lists contained in this class as a fall-back when likelihoods are not available.
 * This comes up in the {@link VariantAnnotator} tool, which annotates variants without the expense of local assembly and pair-HMM.
 *
 */
public class UnfilledReadsLikelihoods<A extends Allele> extends AlleleLikelihoods<GATKRead, A> {

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
            pileupMap.put(samples.getSample(i), ReadPileup.locToReadsPileup(sampleEvidence(i), loc));
        }
        stratifiedPileups = Collections.unmodifiableMap(pileupMap);
        return stratifiedPileups;
    }

    // Methods Which Modify Reads that must be turned off
    @Override
    public void changeEvidence(final Map<GATKRead, GATKRead> readRealignments) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }

    @Override
    public void filterPoorlyModeledEvidence(final ToDoubleFunction<GATKRead> log10MinTrueLikelihood) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }

    @Override
    public void addEvidence(final Map<String,List<GATKRead>> readsBySample, final double initialLikelihood) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }

    @Override
    public void contaminationDownsampling(final Map<String, Double> perSampleDownsamplingFraction) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }

    @Override
    public void filterToOnlyOverlappingEvidence(final SimpleInterval location) {
        throw new UnsupportedOperationException("Cannot alter reads in UnfilledReadsLikelihoods object or cached pileups may be rendered inaccurate, please use a normal ReadsLikelihoods object");
    }
}
