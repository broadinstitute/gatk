package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The infinite-population-genotyping assumes that samples belong to individuals taken at random
 * from a very large population that mate at random.
 * <p>
 *     Consequently genotypes calls between samples are totally independent conditional to the frequencies in
 *     the population they coming from. And genotypes should exhibit the ratios expected under HWE.
 * </p>
 * Therefore each sample genotype likelihoods can be considered to
 * be independent from all other samples.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class InfiniteRandomMatingPopulationModel implements GenotypingModel {

    private final int cachePloidyCapacity;
    private final int cacheAlleleCountCapacity;
    private ThreadLocal<GenotypeLikelihoodCalculator[][]> likelihoodCalculators;

    /**
     * Create a new infinite model instance.
     */
    public InfiniteRandomMatingPopulationModel() {
        this(10,50);
    }

    public InfiniteRandomMatingPopulationModel(final int calculatorCachePloidyCapacity, final int calculatorCacheAlleleCapacity) {
        cachePloidyCapacity = calculatorCachePloidyCapacity;
        cacheAlleleCountCapacity = calculatorCachePloidyCapacity;
        likelihoodCalculators = new ThreadLocal<GenotypeLikelihoodCalculator[][]>( ) {
            @Override
            public GenotypeLikelihoodCalculator[][] initialValue() {
                return new GenotypeLikelihoodCalculator[calculatorCachePloidyCapacity][calculatorCacheAlleleCapacity];
            }
        };
    }

    @Override
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles, final GenotypingData<A> data) {
        if (genotypingAlleles == null) {
            throw new IllegalArgumentException("the allele cannot be null");
        }
        if (data == null) {
            throw new IllegalArgumentException("the genotyping data cannot be null");
        }

        final AlleleListPermutation<A> permutation = data.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = AlleleLikelihoodMatrixMapper.newInstance(permutation);

        final int sampleCount = data.numberOfSamples();

        switch (sampleCount) {
            case 0: return noSampleLikelihoods(permutation,data);
            case 1: return singleSampleLikelihoods(genotypingAlleles,data,alleleLikelihoodMatrixMapper);
            default:
                final PloidyModel ploidyModel = data.ploidyModel();
                return ploidyModel.isHomogeneous() ? multiSampleHomogeneousPloidyModelLikelihoods(genotypingAlleles, data, alleleLikelihoodMatrixMapper, sampleCount, ploidyModel)
                    : multiSampleHeterogeneousPloidyModelLikelihoods(genotypingAlleles, data, alleleLikelihoodMatrixMapper, sampleCount, ploidyModel);
        }
    }

    private <A extends Allele> GenotypingLikelihoods<A> noSampleLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                            final GenotypingData<A> data) {
        @SuppressWarnings("unchecked")
        final List<GenotypeLikelihoods> likelihoods = Collections.EMPTY_LIST;
        return new GenotypingLikelihoods<>(genotypingAlleles,data.ploidyModel(), likelihoods);

    }

    private <A extends Allele> GenotypingLikelihoods<A> singleSampleLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                                final GenotypingData<A> data,
                                                                                final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper) {
        final PloidyModel ploidyModel = data.ploidyModel();
        final int samplePloidy = ploidyModel.samplePloidy(0);
        final int alleleCount = genotypingAlleles.numberOfAlleles();
        final GenotypeLikelihoodCalculator likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy,alleleCount);
        final LikelihoodMatrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.apply(data.readLikelihoods().sampleMatrix(0));
        final List<GenotypeLikelihoods> genotypeLikelihoods = Collections.singletonList(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        return new GenotypingLikelihoods<>(genotypingAlleles,ploidyModel,genotypeLikelihoods);
    }

    private GenotypeLikelihoodCalculator getLikelihoodsCalculator(final int samplePloidy, final int alleleCount) {
        if (samplePloidy >= cacheAlleleCountCapacity) {
            return GenotypeLikelihoodCalculators.getInstance(samplePloidy, alleleCount);
        } else if (alleleCount >= cacheAlleleCountCapacity) {
            return GenotypeLikelihoodCalculators.getInstance(samplePloidy, alleleCount);
        }
        final GenotypeLikelihoodCalculator[][] cache = likelihoodCalculators.get();
        final GenotypeLikelihoodCalculator result = cache[samplePloidy][alleleCount];
        return result != null ? result : (cache[samplePloidy][alleleCount] = GenotypeLikelihoodCalculators.getInstance(samplePloidy, alleleCount));
    }

    private <A extends Allele> GenotypingLikelihoods<A> multiSampleHeterogeneousPloidyModelLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                                                       final GenotypingData<A> data,
                                                                                                       final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper,
                                                                                                       final int sampleCount,
                                                                                                       final PloidyModel ploidyModel) {
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.numberOfAlleles();
        for (int i = 0; i < sampleCount; i++) {
            final int samplePloidy = ploidyModel.samplePloidy(i);
            final GenotypeLikelihoodCalculator likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy,alleleCount);
            final LikelihoodMatrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.apply(data.readLikelihoods().sampleMatrix(i));
            genotypeLikelihoods.add(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles,ploidyModel,genotypeLikelihoods);
    }

    private <A extends Allele> GenotypingLikelihoods<A> multiSampleHomogeneousPloidyModelLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                                                     final GenotypingData<A> data,
                                                                                                     final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper,
                                                                                                     final int sampleCount,
                                                                                                     final PloidyModel ploidyModel) {
        final int samplePloidy = ploidyModel.samplePloidy(0);
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.numberOfAlleles();
        final GenotypeLikelihoodCalculator likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy,alleleCount);
        for (int i = 0; i < sampleCount; i++) {
            final LikelihoodMatrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.apply(data.readLikelihoods().sampleMatrix(i));
            genotypeLikelihoods.add(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles,ploidyModel,genotypeLikelihoods);
    }
}
