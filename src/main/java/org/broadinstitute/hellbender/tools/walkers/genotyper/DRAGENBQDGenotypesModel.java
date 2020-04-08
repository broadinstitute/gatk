package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class DRAGENBQDGenotypesModel implements GenotypersModel {
    private static final int DEFAULT_CACHE_PLOIDY_CAPACITY = 10;
    private static final int DEFAULT_CACHE_ALLELE_CAPACITY = 50;

    private final int cacheAlleleCountCapacity;
    private final int cachePloidyCapacity;
    private GenotypeLikelihoodCalculatorDRAGEN[][] likelihoodCalculators;
    private final GenotypeLikelihoodCalculators calculators;
    private final boolean computeBQD;
    private final boolean computeFRD;
    private final int allelePadding;

    // We keep a fallback model in mind... this might want to be adjusted as implementation workds
    private final IndependentSampleGenotypesModel fallbackModel;

    public DRAGENBQDGenotypesModel(final boolean useBQDModel, final boolean useFRDModel, final int allelePadding) { this(DEFAULT_CACHE_PLOIDY_CAPACITY, DEFAULT_CACHE_ALLELE_CAPACITY, useBQDModel, useFRDModel, allelePadding); }

    /**
     *  Initialize model with given maximum allele count and ploidy for caching
     */
    public DRAGENBQDGenotypesModel(final int calculatorCachePloidyCapacity, final int calculatorCacheAlleleCapacity,
                                   final boolean useBQDModel, final boolean useFRDModel, final int allelePadding) {
        fallbackModel = new IndependentSampleGenotypesModel(calculatorCachePloidyCapacity, calculatorCacheAlleleCapacity);
        cachePloidyCapacity = calculatorCachePloidyCapacity;
        cacheAlleleCountCapacity = calculatorCacheAlleleCapacity;
        likelihoodCalculators = new GenotypeLikelihoodCalculatorDRAGEN[calculatorCachePloidyCapacity][calculatorCacheAlleleCapacity];
        calculators = new GenotypeLikelihoodCalculators();
        this.computeBQD = useBQDModel;
        this.computeFRD = useFRDModel;
        this.allelePadding = allelePadding;
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
        final int variantOffset = data.readLikelihoods().getSubsettedGenomicLoc().getStart() + allelePadding;

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
            List<DragenReadContainer> strandForward = new ArrayList<>();
            List<DragenReadContainer>  strandReverse = new ArrayList<>();

            ////TODO BIG GIANT TODO, THIS IS WRONG!!!! READS WITH INDELS ARE GOING TO BE SORTED INCORRECTLY HERE!!!!!!!!!!!! NEED TO ROLL MY OWN CLIPPING MANAGING CODE.......
            for (int j = 0; j < readsForSample.size(); j++) {
                final GATKRead readForSample = ReadClipper.revertSoftClippedBases(readsForSample.get(j));
                final Pair<Integer, Boolean> baseOffsetForRead = ReadUtils.getReadCoordinateForReferenceCoordinate(readForSample, variantOffset, true);
                final int indexForSnp = readForSample.getBaseQualityCount() > baseOffsetForRead.getLeft() ?  baseOffsetForRead.getLeft() : -1;

                if (readForSample.isReverseStrand()) {
                    strandReverse.add(new DragenReadContainer(readForSample, indexForSnp, ReadUtils.getStrandedUnclippedStart(readForSample), j));
                } else {
                    strandForward.add(new DragenReadContainer(readForSample, indexForSnp, ReadUtils.getStrandedUnclippedStart(readForSample), j));
                }
            }
            //TODO unsilly this
            //TODO marking this with -1s is silly and probably not the right answer
            for (int j = 0; j < hmmFilteredReadsForSample.size(); j++) {
                final GATKRead filteredReadForSample = ReadClipper.revertSoftClippedBases(hmmFilteredReadsForSample.get(j));
                final Pair<Integer, Boolean> baseOffsetForRead = ReadUtils.getReadCoordinateForReferenceCoordinate(filteredReadForSample, variantOffset, true);
                final int indexForSnp = filteredReadForSample.getBaseQualityCount() > baseOffsetForRead.getLeft() ?  baseOffsetForRead.getLeft() : -1; //This is fixing a horrible off by one bug in the above method, this shoulld be fixed but I don't want to deal with the offtarget effects here

                if (filteredReadForSample.isReverseStrand()) {
                    strandReverse.add(new DragenReadContainer(filteredReadForSample, indexForSnp, ReadUtils.getStrandedUnclippedStart(filteredReadForSample), -1));
                } else {
                    strandForward.add(new DragenReadContainer(filteredReadForSample, indexForSnp, ReadUtils.getStrandedUnclippedStart(filteredReadForSample), -1));
                }
            }
            strandForward.sort(new ReadFeatherEndForwardComparitor());
            strandReverse.sort(new ReadFeatherEndRevereseComparitor());

            // Compute default liklihoods as normal (before we go ahead and alter the liklihoods for the call)
            final int samplePloidy = ploidyModel.samplePloidy(sampleIndex);

            // get a new likelihoodsCalculator if this sample's ploidy differs from the previous sample's
            if (samplePloidy != likelihoodsCalculator.ploidy()) {
                likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy, alleleCount);
            }

            // this is the data array for the read liklihoods without any trouble
            final LikelihoodMatrix<GATKRead, A> sampleLikelihoods = alleleLikelihoodMatrixMapper.mapAlleles(data.readLikelihoods().sampleMatrix(sampleIndex));
            final double[] ployidyModelGenotypeLikelihoods = likelihoodsCalculator.rawGenotypeLikelihoods(sampleLikelihoods);

            System.out.println("\n Vanilla resutls:");
            System.out.println(Arrays.toString(ployidyModelGenotypeLikelihoods));

            // TODO these must be instantiated as something real
            double[] BQDCallResults = null;
            double[] FRDCallResults = null;

            if (computeBQD) { // TODO this will become a switch to do frd work or bqd work calling out to the things
                double forwardHomopolymerAdjustment = FRDBQDUtils.computeForwardHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent);
                double reverseHomopolymerAdjustment = FRDBQDUtils.computeReverseHomopolymerAdjustment(paddedReference, offsetForRefIntoEvent);
                BQDCallResults = likelihoodsCalculator.calculateBQDLikelihoods(sampleLikelihoods, strandForward, strandReverse, forwardHomopolymerAdjustment, reverseHomopolymerAdjustment, calculators);
                System.out.println("BQD results:");
                System.out.println(Arrays.toString(BQDCallResults));
            }
            if (computeFRD) { // TODO this will become a switch to do frd work or bqd work calling out to the things
                FRDCallResults = likelihoodsCalculator.calculateFRDLikelihoods(sampleLikelihoods,
                        Stream.of(strandForward, strandReverse).flatMap(list -> list.stream()).filter(r -> r.getIndexInLikelihoodsObject() != -1).collect(Collectors.toList()), // We filter out the HMM filtered reads as they do not apply to FRD
                        34.77, calculators);
                System.out.println("FRD results:");
                System.out.println(Arrays.toString(FRDCallResults));
            }

            //make synthesized likelihoods object (NOTE that we can do this since for invalid model GT fields we simply infinity out the result in the array)
            for (int gt = 0; gt < ployidyModelGenotypeLikelihoods.length; gt++) {
                if (computeBQD) {
                    ployidyModelGenotypeLikelihoods[gt] = Math.max(ployidyModelGenotypeLikelihoods[gt], BQDCallResults[gt]);
                }
                if (computeFRD) {
                    ployidyModelGenotypeLikelihoods[gt] = Math.max(ployidyModelGenotypeLikelihoods[gt], FRDCallResults[gt]);
                }
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

    /**
     * This helper class is used to store the necessary data in order to sort a read based on its BQD "feather end"
     */
    static class DragenReadContainer {
        final public GATKRead underlyingRead;
        final int offsetIntoReadForBaseQuality;
        final int unclippedEnd;
        final int indexInLikelihoodsObject;

        // Transient value used to store thresholds for FRD
       double phredPFValue = 0;


        public DragenReadContainer(final GATKRead underlyingRead, final int offsetIntoReadForBaseQuality, final int unclippedEnd, final int indexInLikelihoodsObject) {
            this.underlyingRead = underlyingRead;
            this.offsetIntoReadForBaseQuality = offsetIntoReadForBaseQuality;
            this.unclippedEnd = unclippedEnd;
            this.indexInLikelihoodsObject = indexInLikelihoodsObject;
        }

        public int getUnclippedPosition() {
            return unclippedEnd;
        }

        public int getIndexInLikelihoodsObject() {
            return indexInLikelihoodsObject;
        }

        public boolean hasValidBaseQuality() {
            return offsetIntoReadForBaseQuality != -1;
        }

        public int getBaseQuality() {
            return underlyingRead.getBaseQuality(offsetIntoReadForBaseQuality);
        }

        public int getMappingQuality() {
            return underlyingRead.getMappingQuality();
        }

        public double getPhredPFValue() {
            return phredPFValue;
        }

        public void setPhredPFValue(double phredPFValue) {
            this.phredPFValue = phredPFValue;
        }

        @Override
        public String toString() {
            return "Read: "+underlyingRead.toString()+" index: "+indexInLikelihoodsObject+" at unclipped end: "+unclippedEnd+" with base quality "+(hasValidBaseQuality() ? getBaseQuality() : -1);
        }

        public boolean isReverseStrand() {
            return underlyingRead.isReverseStrand();
        }
    }

    //MAJOR TODO THIS IS CURRENTLY BASED OFF OF THE REFERENCE UNCLIPPED START AND NOT THE BASES IN THE READ CONSEQUENTLY AT SITES WITH
    //      TODO INDELS PRESENT WE ARE GOING TO BE LOOKING AT THE WRONG OFFSETS FOR THIS SORT... a minor issue but still...

    // Orders the reads based on the number of bases there are to the left of the fatherEndComparisonLocation as aligned according to the cigar
    // NOTE: here we compare the un-hardclipped edges for these reads as the model itself cares about the cycle count of the sequencer, and
    //       importantly this saves us having the thread the original alignment of these reads to this level, since by this point we have trimmed
    //       the reads twice, once to the active region with padding and again to the callable region within the active window and in both of these
    //       cases we have deleted bases with hardclips.
    public class ReadFeatherEndForwardComparitor implements Comparator<DragenReadContainer>, Serializable {
        private static final long serialVersionUID = 1L;
        /**
         * Evaluate first the number of bases to the end of the read and follow that by the base quality
         */
        @Override
        public int compare(final DragenReadContainer read1, final DragenReadContainer read2) {
            //NOTE: here we want the reads to wind up in ascending order by unclipped position because the unclipped position should be on the left
            int diffVal = read1.getUnclippedPosition() - read2.getUnclippedPosition();
            if (diffVal == 0) {
                diffVal = (read1.hasValidBaseQuality() ? read1.getBaseQuality() : 0)
                        - (read2.hasValidBaseQuality() ? read2.getBaseQuality() : 0);
            }
            return diffVal;
        }
    }

    // Orders the reads based on the number of bases in the read that occur before the fatherEndComparisonLocation as aligned according to the cigar
    // NOTE: here we compare the un-hardclipped edges for these reads as the model itself cares about the cycle count of the sequencer, and
    //       importantly this saves us having the thread the original alignment of these reads to this level, since by this point we have trimmed
    //       the reads twice, once to the active region with padding and again to the callable region within the active window and in both of these
    //       cases we have deleted bases with hardclips.
    public class ReadFeatherEndRevereseComparitor implements Comparator<DragenReadContainer>, Serializable {
        private static final long serialVersionUID = 1L;
        /**
         * Evaluate first the number of bases to the end of the read and follow that by the base quality
         */
        @Override
        public int compare(final DragenReadContainer read1, final DragenReadContainer read2) {
            //NOTE: here we want the reads to wind up in decending order by unclipped position because the unclipped position should be on the left
            int diffVal = read2.getUnclippedPosition() - read1.getUnclippedPosition();
            if (diffVal==0) {
                //TODO verify this lines up with the sort in DRAGBQD
                diffVal = (read1.hasValidBaseQuality() ? read1.getBaseQuality() : 0)
                        - (read2.hasValidBaseQuality() ? read2.getBaseQuality() : 0);
            }
            return diffVal;
        }

    }
}
