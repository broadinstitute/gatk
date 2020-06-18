package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.transformers.DRAGENMappingQualityReadTransformer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrParams;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrReferenceSTRs;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This is the DRAGEN-GATK genotyper model. This class manages the logic for likelihoods calculation between the IndependentSamplesGenotyperModel
 * and the two DRAGEN genotypes models (FRD and BQD).
 *
 * In order to access the FRD and BQD genotypes model simply initialize this class with the settings enabled, then call
 * {@link #calculateLikelihoods} to return a the modified likelihoods array with the relevant models applied. FRD and BQD
 * both work by computing alternate likelihoods scores for the array array corresponding to homozygous allele combinations
 * and choosing the best score for each element of the likelihoods array from among the different models. This has the net
 * effect of penalizing heterozygous allele combinations since they do not have their likelihoods modified.
 *
 */
public class DRAGENGenotypesModel implements GenotypersModel {
    private static final int DEFAULT_CACHE_PLOIDY_CAPACITY = 10;
    private static final int DEFAULT_CACHE_ALLELE_CAPACITY = 50;
    // Flat SNP het prior to use for genotyping
    public static final double FLAT_SNP_HET_PRIOR = 34.77;

    private final int cacheAlleleCountCapacity;
    private final int cachePloidyCapacity;
    private GenotypeLikelihoodCalculatorDRAGEN[][] likelihoodCalculators;
    private final GenotypeLikelihoodCalculators calculators;
    private final boolean computeBQD;
    private final boolean computeFRD;
    private final int allelePadding;
    private final int maxEffectiveDepthAdjustment;
    private final DragstrParams dragstrParams;

    //Debug stream to be managed by the HaplotypeCallerEngine
    private PrintStream genotyperDebugStream = null;

    public DRAGENGenotypesModel(final boolean useBQDModel, final boolean useFRDModel, final int allelePadding, final int maxEffectiveDepthAdjustment, final DragstrParams dragstrParams) { this(DEFAULT_CACHE_PLOIDY_CAPACITY, DEFAULT_CACHE_ALLELE_CAPACITY, useBQDModel, useFRDModel, allelePadding, maxEffectiveDepthAdjustment,  dragstrParams); }

    /*
     *  Initialize model with given maximum allele count and ploidy for caching
     */
    public DRAGENGenotypesModel(final int calculatorCachePloidyCapacity, final int calculatorCacheAlleleCapacity,
                                final boolean useBQDModel, final boolean useFRDModel, final int allelePadding,
                                final int maxEffectiveDepthAdjustment, final DragstrParams dragstrParams) {
        cachePloidyCapacity = calculatorCachePloidyCapacity;
        cacheAlleleCountCapacity = calculatorCacheAlleleCapacity;
        likelihoodCalculators = new GenotypeLikelihoodCalculatorDRAGEN[calculatorCachePloidyCapacity][calculatorCacheAlleleCapacity];
        calculators = new GenotypeLikelihoodCalculators();
        this.computeBQD = useBQDModel;
        this.computeFRD = useFRDModel;
        this.allelePadding = allelePadding;
        this.maxEffectiveDepthAdjustment = maxEffectiveDepthAdjustment;
        this.dragstrParams = dragstrParams;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles, final GenotypingData<A> data, byte[] paddedReference, int offsetForRefIntoEvent, final DragstrReferenceSTRs dragstrs) {
        Utils.nonNull(genotypingAlleles, "the allele cannot be null");
        Utils.nonNull(data, "the genotyping data cannot be null");

        //Get the prior for the
        double api;
        if (dragstrs !=  null) {
            final int period = dragstrs.period(offsetForRefIntoEvent + 1 );
            final int repeats = dragstrs.repeatLength(offsetForRefIntoEvent + 1);
            api = dragstrParams.api(period, repeats);
            if (genotyperDebugStream != null) {
                genotyperDebugStream.println("API found: " + api + " with period used: " + period + "  and repeats: " + repeats);
            }
        } else {
            if (genotyperDebugStream != null) {
                genotyperDebugStream.println("No API from DRAGStrs found, falling back on snp het prior for indels");
            }
            api = FLAT_SNP_HET_PRIOR;
        }


        final AlleleListPermutation<A> permutation = data.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = new AlleleLikelihoodMatrixMapper(permutation);

        final int sampleCount = data.numberOfSamples();
        final PloidyModel ploidyModel = data.ploidyModel();
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.numberOfAlleles();
        final int variantOffset = data.readLikelihoods().getVariantCallingSubsetApplied().getStart() + allelePadding;

        GenotypeLikelihoodCalculatorDRAGEN likelihoodsCalculator = getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount); //TODO this needs to change
        for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {

            ///////////////////////////////////////////////////////////////////////////
            ///// PREPROCESSING FOR BQD
            ///////////////////////////////////////////////////////////////////////////

            // Separating the reads by their strand and sorting them appropriately.
            List<GATKRead> readsForSample = data.readLikelihoods().sampleEvidence(sampleIndex);
            List<GATKRead> hmmFilteredReadsForSample = data.readLikelihoods().filteredSampleEvidence(sampleIndex);
            // These objects are intended to store 3 things, the read, the inner (middle) int stores the offset into the read of the base in question, and the outer int stores the index of the read per sample
            List<DragenReadContainer> strandForward = new ArrayList<>();
            List<DragenReadContainer>  strandReverse = new ArrayList<>();

            ////TODO reads with indels preceding the variant in question might have their cycle counts mismatched, its unclear whether dragen handles this case based on debug outputs
            for (int j = 0; j < readsForSample.size(); j++) {
                final GATKRead readForSample = readsForSample.get(j);
                final int indexForSnp = ReadUtils.getReadIndexForReferenceCoordinate(readForSample, variantOffset).getLeft();

                if (readForSample.isReverseStrand()) {
                    strandReverse.add(new DragenReadContainer(readForSample, indexForSnp, ReadUtils.getStrandedUnclippedStart(readForSample), j));
                } else {
                    strandForward.add(new DragenReadContainer(readForSample, indexForSnp, ReadUtils.getStrandedUnclippedStart(readForSample), j));
                }
            }
            for (int j = 0; j < hmmFilteredReadsForSample.size(); j++) {
                final GATKRead filteredReadForSample = hmmFilteredReadsForSample.get(j);
                final int indexForSnp = ReadUtils.getReadIndexForReferenceCoordinate(filteredReadForSample, variantOffset).getLeft();

                if (filteredReadForSample.isReverseStrand()) {
                    strandReverse.add(new DragenReadContainer(filteredReadForSample, indexForSnp, ReadUtils.getStrandedUnclippedStart(filteredReadForSample), -1));
                } else {
                    strandForward.add(new DragenReadContainer(filteredReadForSample, indexForSnp, ReadUtils.getStrandedUnclippedStart(filteredReadForSample), -1));
                }
            }
            strandForward.sort(new ReadFeatherEndForwardComparitor());//data.readLikelihoods().getVariantCallingSubsetApplied().getContig(), variantOffset);
            strandReverse.sort(new ReadFeatherEndRevereseComparitor());

            // Compute default liklihoods as normal (before we go ahead and alter the liklihoods for the call)
            final int samplePloidy = ploidyModel.samplePloidy(sampleIndex);

            // get a new likelihoodsCalculator if this sample's ploidy differs from the previous sample's
            if (samplePloidy != likelihoodsCalculator.ploidy()) {
                likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy, alleleCount);
            }
            if (genotyperDebugStream != null) {
                likelihoodsCalculator.addGenotyperDebugOutputStream(genotyperDebugStream);
            }

            // this is the data array for the read liklihoods without any trouble
            final LikelihoodMatrix<GATKRead, A> sampleLikelihoods = alleleLikelihoodMatrixMapper.mapAlleles(data.readLikelihoods().sampleMatrix(sampleIndex));
            final double[] ployidyModelGenotypeLikelihoods = likelihoodsCalculator.rawGenotypeLikelihoods(sampleLikelihoods);

            if (genotyperDebugStream != null) {
                genotyperDebugStream.println("\n Standard Genotyping Resutls:");
                genotyperDebugStream.println(Arrays.toString(ployidyModelGenotypeLikelihoods));
            }

            double[] BQDCallResults = null;
            double[] FRDCallResults = null;
            if (computeBQD) {
                BQDCallResults = likelihoodsCalculator.calculateBQDLikelihoods(sampleLikelihoods, strandForward, strandReverse, paddedReference, offsetForRefIntoEvent, calculators);
                if (genotyperDebugStream != null) {
                    genotyperDebugStream.println("BQD results:");
                    genotyperDebugStream.println(Arrays.toString(BQDCallResults));
                }
            }
            if (computeFRD) {
                FRDCallResults = likelihoodsCalculator.calculateFRDLikelihoods(sampleLikelihoods, ployidyModelGenotypeLikelihoods,
                        Stream.of(strandForward, strandReverse).flatMap(Collection::stream).collect(Collectors.toList()), // We filter out the HMM filtered reads as they do not apply to FRD
                        FLAT_SNP_HET_PRIOR, api, maxEffectiveDepthAdjustment, calculators);
                if (genotyperDebugStream != null) {
                    genotyperDebugStream.println("FRD results:");
                    genotyperDebugStream.println(Arrays.toString(FRDCallResults));
                }
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
            if (genotyperDebugStream != null) {
                genotyperDebugStream.println("merged matrix:");
                genotyperDebugStream.println(Arrays.toString(ployidyModelGenotypeLikelihoods));
            }
        }
        return new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods);
    }

    // Adding the debug stream managed by the haplotype caller engine
    @Override
    public void addDebugOutStream(final PrintStream debugStream) {
        this.genotyperDebugStream = debugStream;
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

    /**
     * This helper class is used to store the necessary data in order to sort a read based on its BQD "feather end" as
     * well as information relevant to re-associate the read with its position in the AlleleLikelihoods object arrays.
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

        public boolean wasFilteredByHMM() {
            return this.getIndexInLikelihoodsObject() == -1;
        }

        public boolean hasValidBaseQuality() {
            return offsetIntoReadForBaseQuality != -1;
        }

        public int getBaseQuality() {
            return underlyingRead.getBaseQuality(offsetIntoReadForBaseQuality);
        }

        public int getForwardsFeatherEnd() {
            return (underlyingRead.getSoftStart() - underlyingRead.getUnclippedStart()) + offsetIntoReadForBaseQuality;
        }

        public int getReverseFeatherEnd() {
            return (underlyingRead.getUnclippedEnd() - underlyingRead.getSoftEnd()) + (underlyingRead.getLength() - offsetIntoReadForBaseQuality);
        }

        public double getPhredScaledMappingQuality() {
            return DRAGENMappingQualityReadTransformer.mapValue(underlyingRead.getMappingQuality());
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

            int diffVal =  read2.getForwardsFeatherEnd() - read1.getForwardsFeatherEnd();
//
//            int diffVal = read1.getUnclippedPosition() - read2.getUnclippedPosition();
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
            int diffVal = read2.getReverseFeatherEnd() - read1.getReverseFeatherEnd();
//            int diffVal = read2.getUnclippedPosition() - read1.getUnclippedPosition();
            if (diffVal==0) {
                diffVal = (read1.hasValidBaseQuality() ? read1.getBaseQuality() : 0)
                        - (read2.hasValidBaseQuality() ? read2.getBaseQuality() : 0);
            }
            return diffVal;
        }

    }
}
