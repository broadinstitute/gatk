package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

/**
 * Encapsulates the data use to make the genotype calls.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GenotypingData<A extends Allele> implements SampleList, AlleleList<A> {

    private final PloidyModel ploidyModel;

    private final ReadLikelihoods<A> likelihoods;

    /**
     * Constructs a new genotyping-data collection providing the ploidy model to apply to the input model
     * and the read-likelihoods collection.
     *
     *
     * @param ploidyModel the ploidy model.
     * @param likelihoods the read-likelihoods collection.
     *
     * @throws IllegalArgumentException if either {@code ploidyModel} or {@code likelihoods} is {@code null},
     *   or they are not compatible in terms of the samples they contain; their lists must match.
     */
    public GenotypingData(final PloidyModel ploidyModel, final ReadLikelihoods<A> likelihoods) {
        if (ploidyModel == null) {
            throw new IllegalArgumentException("the ploidy model cannot be null");
        }
        if (likelihoods == null) {
            throw new IllegalArgumentException("the likelihood object cannot be null");
        }
        this.ploidyModel = ploidyModel;
        this.likelihoods = likelihoods;
        if (!ploidyModel.asListOfSamples().equals(likelihoods.asListOfSamples())) {
            throw new IllegalArgumentException("sample list are different between ploidy-model and read-likelihood collection, perhaps just the order");
        }
    }

    /**
     * Returns the ploidy model that corresponds to the data provided.
     * @return never {@code null}.
     */
    public PloidyModel ploidyModel() {
        return ploidyModel;
    }

    @Override
    public int numberOfSamples() {
        return ploidyModel.numberOfSamples();
    }

    @Override
    public int indexOfSample(final String sample) {
        return ploidyModel.indexOfSample(sample);
    }

    @Override
    public String getSample(final int sampleIndex) {
        return ploidyModel.getSample(sampleIndex);
    }

    /**
     * Returns read-likelihoods to use for genotyping.
     * @return never {@code null}.
     */
    public ReadLikelihoods<A> readLikelihoods() {
        return likelihoods;
    }

    @Override
    public int numberOfAlleles() {
        return likelihoods.numberOfAlleles();
    }

    @Override
    public int indexOfAllele(final A allele) {
        return likelihoods.indexOfAllele(allele);
    }

    @Override
    public A getAllele(final int index) {
        return likelihoods.getAllele(index);
    }
}
