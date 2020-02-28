package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DRAGENBQDGenotypesModel implements GenotypersModel {
    private static final int DEFAULT_CACHE_PLOIDY_CAPACITY = 10;
    private static final int DEFAULT_CACHE_ALLELE_CAPACITY = 50;

    private final int cacheAlleleCountCapacity;
    private final int cachePloidyCapacity;
    private GenotypeLikelihoodCalculatorDRAGEN[][] likelihoodCalculators;
    private final GenotypeLikelihoodCalculators calculators;
    private final boolean computeBQD;
    private final boolean computeFRD;

    // We keep a fallback model in mind... this might want to be adjusted as implementation workds
    private final IndependentSampleGenotypesModel fallbackModel;

    public DRAGENBQDGenotypesModel(final boolean useBQDModel, final boolean useFRDModel) { this(DEFAULT_CACHE_PLOIDY_CAPACITY, DEFAULT_CACHE_ALLELE_CAPACITY, useBQDModel, useFRDModel); }

    /**
     *  Initialize model with given maximum allele count and ploidy for caching
     */
    public DRAGENBQDGenotypesModel(final int calculatorCachePloidyCapacity, final int calculatorCacheAlleleCapacity,
                                   final boolean useBQDModel, final boolean useFRDModel) {
        fallbackModel = new IndependentSampleGenotypesModel(calculatorCachePloidyCapacity, calculatorCacheAlleleCapacity);
        cachePloidyCapacity = calculatorCachePloidyCapacity;
        cacheAlleleCountCapacity = calculatorCacheAlleleCapacity;
        likelihoodCalculators = new GenotypeLikelihoodCalculatorDRAGEN[calculatorCachePloidyCapacity][calculatorCacheAlleleCapacity];
        calculators = new GenotypeLikelihoodCalculators();
        this.computeBQD = useBQDModel;
        this.computeFRD = useFRDModel;
    }

//    @Override TODO unify the two calling infrastructures
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles, final GenotypingData<A> data, byte[] paddedReference, int offsetForRefIntoEvent) {
        Utils.nonNull(genotypingAlleles, "the allele cannot be null");
        Utils.nonNull(data, "the genotyping data cannot be null");

        // for right now, don't handle any deletions whatsoever
        // Also for right now lets not worry too much abotu alleles.
        if (FRDBQDUtils.containsInsertionOrDeletion(genotypingAlleles) || data.numberOfAlleles() > 3) {
            return fallbackModel.calculateLikelihoods(genotypingAlleles, data);
        }

        final AlleleListPermutation<A> permutation = data.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = new AlleleLikelihoodMatrixMapper(permutation);

        final int sampleCount = data.numberOfSamples();
        final PloidyModel ploidyModel = data.ploidyModel();
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.numberOfAlleles();
        final int variantOffset = data.readLikelihoods().getSubsettedGenomicLoc().getStart() + HaplotypeCallerGenotypingEngine.ALLELE_EXTENSION;


//        GenotypeLikelihoodCalculator likelihoodsCalculator = sampleCount > 0 ? getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount) : null;
        GenotypeLikelihoodCalculatorDRAGEN likelihoodsCalculator = getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount); //TODO this needs to change
        for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {

            ///////////////////////////////////////////////////////////////////////////
            ///// PREPROCESSING FOR BQD
            ///////////////////////////////////////////////////////////////////////////
            //TODO this will need ohmptimation, this is probably too much object creation for in here

            // Separating the reads by their strand and sorting them appropriately.
            List<GATKRead> readsForSample = data.readLikelihoods().sampleEvidence(sampleIndex);
            List<GATKRead> hmmFilteredReadsForSample = data.readLikelihoods().filteredSampleEvidence(sampleIndex);
            // These objects are intended to store 3 things, the read, the inner (middle) int stores the offset into the read of the base in question, and the outer int stores the index of the read per sample
            List<Pair<Pair<GATKRead,Integer>,Integer>> strandForward = new ArrayList<>();
            List<Pair<Pair<GATKRead,Integer>,Integer>>  strandReverse = new ArrayList<>();
            for (int j = 0; j < readsForSample.size(); j++) {
                // TODO figure out what to do with overlapping deletions....
                final Pair<Integer, Boolean> baseOffsetForRead = ReadUtils.getReadCoordinateForReferenceCoordinate(readsForSample.get(j), variantOffset, true);
                if (readsForSample.get(j).isReverseStrand()) {
                    strandReverse.add(Pair.of(
                            Pair.of(readsForSample.get(j),baseOffsetForRead.getLeft()), j));
                } else {
                    strandForward.add(Pair.of(
                            Pair.of(readsForSample.get(j),baseOffsetForRead.getLeft()), j));
                }
            }
            //TODO unsilly this
            //TODO marking this with -1s is silly and probably not the right answer
            for (int j = 0; j < hmmFilteredReadsForSample.size(); j++) {
                final Pair<Integer, Boolean> baseOffsetForRead = ReadUtils.getReadCoordinateForReferenceCoordinate(hmmFilteredReadsForSample.get(j), variantOffset, true);
                if (hmmFilteredReadsForSample.get(j).isReverseStrand()) {
                    strandReverse.add(Pair.of(
                            Pair.of(hmmFilteredReadsForSample.get(j),baseOffsetForRead.getLeft()), -1));
                } else {
                    strandForward.add(Pair.of(
                            Pair.of(hmmFilteredReadsForSample.get(j),baseOffsetForRead.getLeft()), -1));
                }
            }
            strandForward.sort(new FRDBQDUtils.ReadFeatherEndForwardComparitor());
            strandReverse.sort(new FRDBQDUtils.ReadFeatherEndRevereseComparitor());

            // Compute default liklihoods as normal (before we go ahead and alter the liklihoods for the call)
            final int samplePloidy = ploidyModel.samplePloidy(sampleIndex);

            // get a new likelihoodsCalculator if this sample's ploidy differs from the previous sample's
            if (samplePloidy != likelihoodsCalculator.ploidy()) {
                likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy, alleleCount);
            }

            // this is the data array for the read liklihoods without any trouble
            final LikelihoodMatrix<GATKRead, A> sampleLikelihoods = alleleLikelihoodMatrixMapper.mapAlleles(data.readLikelihoods().sampleMatrix(sampleIndex));
            final double[] ployidyModelGenotypeLikelihoods = likelihoodsCalculator.rawGenotypeLikelihoods(sampleLikelihoods);

            System.out.println("Genotyping model results for alleles before being modified");
            System.out.println(ployidyModelGenotypeLikelihoods.toString());

            // TODO these must be instantiated as something real
            double[] BQDCallResults = null;
            double[] FRDCallResults = null;

            if (computeBQD) { // TODO this will become a switch to do frd work or bqd work calling out to the things
                double forwardHomopolymerAdjustment = FRDBQDUtils.computeForwardHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent);
                double reverseHomopolymerAdjustment = FRDBQDUtils.computeReverseHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent);
                BQDCallResults = likelihoodsCalculator.calculateBQDLikelihoods(sampleLikelihoods, strandForward, strandReverse, forwardHomopolymerAdjustment, reverseHomopolymerAdjustment, calculators);
            }

            System.out.println("Genotyping model results for genotypes given BQD results");
            System.out.println(Arrays.asList(BQDCallResults));

            if (computeFRD) { // TODO this will become a switch to do frd work or bqd work calling out to the things
//                FRDCallResults = likelihoodsCalculator.calculateFRDLikelihoods(sampleLikelihoods, strandForward, strandReverse, forwardHomopolymerAdjustment, reverseHomopolymerAdjustment);
            }

            System.out.println("Genotyping model results for genotypes given FRD results");
            System.out.println(Arrays.asList(FRDCallResults));

            //make synthesized likelihoods object (NOTE that we can do this since for invalid model GT fields we simply infinity out the result in the array)
            for (int gt = 0; gt < ployidyModelGenotypeLikelihoods.length; sampleIndex++) {
                ployidyModelGenotypeLikelihoods[gt] = Math.min(ployidyModelGenotypeLikelihoods[gt], Math.min(BQDCallResults[gt], FRDCallResults[gt]));
            }

            // this is what the work actually is, after we have computed a few things
            genotypeLikelihoods.add(GenotypeLikelihoods.fromLog10Likelihoods(ployidyModelGenotypeLikelihoods));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods);
    }

    private GenotypeLikelihoodCalculatorDRAGEN getLikelihoodsCalculator(final int samplePloidy, final int alleleCount) {
        if (samplePloidy >= cachePloidyCapacity || alleleCount >= cacheAlleleCountCapacity) {
            return calculators.getInstanceDRAGEN(samplePloidy, alleleCount);
        }
        final GenotypeLikelihoodCalculatorDRAGEN result = likelihoodCalculators[samplePloidy][alleleCount];
        if (result != null) {
            return result;
        } else {
            final GenotypeLikelihoodCalculatorDRAGEN newOne = calculators.getInstanceDRAGEN(samplePloidy, alleleCount);
            likelihoodCalculators[samplePloidy][alleleCount] = newOne;
            return newOne;
        }
    }

    // TODO this can be done away with at some future date
    @Override
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(AlleleList<A> genotypingAlleles, GenotypingData<A> data) {
        return null;
    }
}
