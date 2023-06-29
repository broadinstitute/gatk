package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Refines variant breakpoints using split read evidence.
 *
 * The start and end of the breakpoint are tested independently. At each end we perform a series of Poisson
 * tests using the following model:
 *
 *  Given:
 *      d_i : depth of sample i
 *      c_i : split read count of sample i
 *
 *   Define a Poisson model of split read counts:
 *      r_i = c_i / d_i : normalized split read count of sample i
 *      m_c : median carrier r_i
 *      m_b : median background (non-carrier) r_i
 *      mu : mean depth of all samples
 *
 *      lambda = mu * m_c : expected carrier count
 *      X ~ Poisson(lambda) : carrier count model
 *      x_b = round(mu * m_b) : adjusted median background count
 *
 *   Calculate probability of observing the background count:
 *      p = P(X < x_b)
 *
 *   We then select the site with the lowest score. Breakpoint end positions are restricted by a lowerbound that depends
 *   on the refined start position (see {@link BreakpointRefiner#getEndLowerBound(SVCallRecord, int)}.
 */
public class BreakpointRefiner {

    private final Map<String,Double> sampleCoverageMap;
    private final SAMSequenceDictionary dictionary;
    /**
     * Number bases that split read positions can pass by the original breakpoint, when left-clipped
     * reads are to the left of the breakpoint and right-clipped reads to the right. Applies only to INS/DEL.
     */
    protected int maxSplitReadCrossDistance;
    protected int representativeDepth;

    public static final int DEFAULT_MAX_CROSS_DISTANCE = 200;
    public static final int MAX_QUAL = 99;

    /**
     * @param sampleCoverageMap map with (sample id, per-base sample coverage) entries
     * @param dictionary reference dictionary
     */
    public BreakpointRefiner(final Map<String, Double> sampleCoverageMap, int maxSplitReadCrossDistance,
                             final SAMSequenceDictionary dictionary) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.dictionary = Utils.nonNull(dictionary);
        this.maxSplitReadCrossDistance = maxSplitReadCrossDistance;
        this.representativeDepth = EvidenceStatUtils.computeRepresentativeDepth(sampleCoverageMap.values());
    }

    /**
     * Performs refinement on one side of a breakpoint
     *
     * @param sortedEvidence split read evidence to test, sorted by position
     * @param strand split read evidence strand
     * @param carrierSamples carrier sample ids
     * @param backgroundSamples background sample ids
     * @param representativeDepth normalization depth
     * @param defaultContig contig to use if test cannot be performed (no evidence or carriers)
     * @param defaultPosition position to use if test cannot be performed (no evidence or carriers)
     * @return pair containing site with refined breakpoints and probability (null if no evidence or carriers)
     */
    protected static SplitReadSite refineSplitReadSite(final List<SplitReadEvidence> sortedEvidence,
                                                       final boolean strand,
                                                       final Collection<String> carrierSamples,
                                                       final Collection<String> backgroundSamples,
                                                       final Map<String, Double> sampleCoverageMap,
                                                       final int representativeDepth,
                                                       final String defaultContig,
                                                       final int defaultPosition) {
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(carrierSamples),
                "One or more carrier samples not found in sample coverage map");
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(backgroundSamples),
                "One or more non-carrier samples not found in sample coverage map");

        // Default case
        if (sortedEvidence.isEmpty() || carrierSamples.isEmpty()) {
            return new SplitReadSite(defaultContig, defaultPosition, strand, Collections.emptyMap(), null);
        }

        EvidenceStatUtils.PoissonTestResult minPResult = null;
        Integer minDistance = null;
        Integer minPPosition = null;
        String minPContig = null;
        Map<String, Integer> minPSampleCounts = null;
        int position;
        Map<String, Integer> sampleCounts = new HashMap<>();
        for (int i = 0; i < sortedEvidence.size(); i++) {
            final SplitReadEvidence e = sortedEvidence.get(i);
            position = e.getStart();
            sampleCounts.put(e.getSample(), e.getCount());
            if (i == sortedEvidence.size() - 1 || sortedEvidence.get(i + 1).getStart() != position) {
                final EvidenceStatUtils.PoissonTestResult result = EvidenceStatUtils.calculateOneSamplePoissonTest(
                        sampleCounts, carrierSamples, backgroundSamples, sampleCoverageMap, representativeDepth
                );
                final int dist = Math.abs(position - defaultPosition);
                if (minPResult == null || result.getP() < minPResult.getP() || (result.getP() == minPResult.getP() && dist < minDistance)) {
                    minPResult = result;
                    minPPosition = position;
                    minPContig = e.getContig();
                    minPSampleCounts = sampleCounts;
                    minDistance = dist;
                }
                sampleCounts = new HashMap<>();
            }
        }
        return new SplitReadSite(minPContig, minPPosition, strand, minPSampleCounts, minPResult);
    }

    /**
     * Performs breakend refinement for a call
     *
     * @param record with split read evidence
     * @return record with new breakpoints
     */
    public RefineResult testRecord(final SVCallRecord record,
                                   final List<SplitReadEvidence> startEvidence,
                                   final List<SplitReadEvidence> endEvidence,
                                   final Set<String> carrierSamples,
                                   final Set<String> backgroundSamples,
                                   final DiscordantPairEvidenceTester.DiscordantPairTestResult discordantPairResult) {
        Utils.nonNull(record);
        SVCallRecordUtils.validateCoordinatesWithDictionary(record, dictionary);
        Utils.validateArg(record.getStrandA() != null, "Record has null strand A");
        Utils.validateArg(record.getStrandB() != null, "Record has null strand B");
        SplitReadSite refinedFirstSite;
        SplitReadSite refinedSecondSite;
        if (!record.isIntrachromosomal()) {
            // Interchromosomal variants, just refine without any checks
            refinedFirstSite = refineSplitReadSite(startEvidence, record.getStrandA(), carrierSamples,
                    backgroundSamples, sampleCoverageMap, representativeDepth, record.getContigA(), record.getPositionA());
            refinedSecondSite = refineSplitReadSite(endEvidence, record.getStrandB(), carrierSamples,
                    backgroundSamples, sampleCoverageMap, representativeDepth, record.getContigB(), record.getPositionB());
        } else if (record.getStrandA() == record.getStrandB()) {
            // Case of intrachromosomal and matching strands, need to ensure start/end don't end up the same
            refinedFirstSite = refineSplitReadSite(startEvidence, record.getStrandA(), carrierSamples,
                    backgroundSamples, sampleCoverageMap, representativeDepth, record.getContigA(), record.getPositionA());
            // Filter out evidence at the refined start position so that we don't double-test
            final int refinedStartPosition = refinedFirstSite.getPosition();
            final List<SplitReadEvidence> validEndEvidence = endEvidence.stream()
                    .filter(e -> e.getStart() != refinedStartPosition).collect(Collectors.toList());
            refinedSecondSite = refineSplitReadSite(validEndEvidence, record.getStrandB(), carrierSamples,
                    backgroundSamples, sampleCoverageMap, representativeDepth, record.getContigB(), record.getPositionB());
        } else {
            // Intrachromosomal and non-matching strands
            refinedFirstSite = refineSplitReadSite(startEvidence, record.getStrandA(), carrierSamples,
                    backgroundSamples, sampleCoverageMap, representativeDepth, record.getContigA(), record.getPositionA());
            refinedSecondSite = refineSplitReadSite(endEvidence, record.getStrandB(), carrierSamples,
                    backgroundSamples, sampleCoverageMap, representativeDepth, record.getContigB(), record.getPositionB());

            // Check if refined coordinates are valid. If not, choose the better one and refine the other again
            final int endLowerBound = getEndLowerBound(record, refinedFirstSite.getPosition());
            if (refinedSecondSite.getPosition() < endLowerBound) {
                if (refinedFirstSite.getP() < refinedSecondSite.getP()) {
                    // Start site had more significant result, so recompute end site with valid boundaries
                    final List<SplitReadEvidence> validEndEvidence = filterSplitReadSitesLowerBound(endEvidence, endLowerBound);
                    final int defaultEndPosition = Math.max(endLowerBound, record.getPositionB());
                    refinedSecondSite = refineSplitReadSite(validEndEvidence, record.getStrandB(), carrierSamples,
                            backgroundSamples, sampleCoverageMap, representativeDepth, record.getContigB(), defaultEndPosition);
                } else {
                    // End site had more significant result, so recompute start site with valid boundaries
                    final int startUpperBound = getStartUpperBound(record, refinedSecondSite.getPosition());
                    final List<SplitReadEvidence> validStartEvidence = filterSplitReadSitesUpperBound(startEvidence, startUpperBound);
                    final int defaulStartPosition = Math.min(startUpperBound, record.getPositionA());
                    refinedFirstSite = refineSplitReadSite(validStartEvidence, record.getStrandA(), carrierSamples,
                            backgroundSamples, sampleCoverageMap, representativeDepth, record.getContigA(), defaulStartPosition);
                }
            }
        }

        // Compute stats on sum of start and end counts
        final EvidenceStatUtils.PoissonTestResult bothsideResult = calculateBothsideTest(refinedFirstSite, refinedSecondSite,
                carrierSamples, backgroundSamples, sampleCoverageMap, representativeDepth);

        EvidenceStatUtils.PoissonTestResult combinedResult = null;
        if (discordantPairResult != null) {
            combinedResult = calculatePESRTest(refinedFirstSite, refinedSecondSite,
                    discordantPairResult, carrierSamples, backgroundSamples,
                    sampleCoverageMap, representativeDepth);
        }

        return new RefineResult(refinedFirstSite, refinedSecondSite, bothsideResult, discordantPairResult, combinedResult);
    }

    public SVCallRecord applyToRecord(final SVCallRecord record,
                                      final RefineResult result) {
        Utils.nonNull(record);
        Utils.nonNull(result);
        final SplitReadSite refinedFirstSite;
        final SplitReadSite refinedSecondSite;
        if (record.isIntrachromosomal() && result.getSecond().getPosition() < result.getFirst().getPosition()) {
            // Swap first and second if out of order
            refinedFirstSite = result.getSecond();
            refinedSecondSite = result.getFirst();
        } else {
            refinedFirstSite = result.getFirst();
            refinedSecondSite = result.getSecond();
        }
        final EvidenceStatUtils.PoissonTestResult bothsideResult = result.getBothsidesResult();

        final int newStart;
        final int newEnd;
        if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
            // By convention, choose the left coordinate, which is guaranteed after the above check
            newStart = refinedFirstSite.getPosition();
            newEnd = newStart;
        } else {
            newStart = refinedFirstSite.getPosition();
            newEnd = refinedSecondSite.getPosition();
        }

        final Integer length = record.getType().equals(GATKSVVCFConstants.StructuralVariantAnnotationType.INS) ? record.getLength() : null;

        final Integer firstQuality = refinedFirstSite.getP() == null || Double.isNaN(refinedFirstSite.getP()) ? null : EvidenceStatUtils.probToQual(refinedFirstSite.getP(), (byte) MAX_QUAL);
        final Integer secondQuality = refinedSecondSite.getP() == null || Double.isNaN(refinedSecondSite.getP()) ? null : EvidenceStatUtils.probToQual(refinedSecondSite.getP(), (byte) MAX_QUAL);
        final Integer totalQuality = Double.isNaN(bothsideResult.getP()) ? null : EvidenceStatUtils.probToQual(bothsideResult.getP(), (byte) MAX_QUAL);
        final Map<String, Object> refinedAttr = new HashMap<>(record.getAttributes());
        refinedAttr.put(GATKSVVCFConstants.FIRST_SPLIT_QUALITY_ATTRIBUTE, firstQuality);
        refinedAttr.put(GATKSVVCFConstants.SECOND_SPLIT_QUALITY_ATTRIBUTE, secondQuality);
        refinedAttr.put(GATKSVVCFConstants.TOTAL_SPLIT_QUALITY_ATTRIBUTE, totalQuality);

        final Integer startCarrierSignal = EvidenceStatUtils.carrierSignalFraction(refinedFirstSite.getCarrierSignal(),
                refinedFirstSite.getBackgroundSignal());
        final Integer endCarrierSignal = EvidenceStatUtils.carrierSignalFraction(refinedSecondSite.getCarrierSignal(),
                refinedSecondSite.getBackgroundSignal());
        final Integer totalCarrierSignal = EvidenceStatUtils.carrierSignalFraction(bothsideResult.getCarrierSignal(),
                bothsideResult.getBackgroundSignal());
        refinedAttr.put(GATKSVVCFConstants.FIRST_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, startCarrierSignal);
        refinedAttr.put(GATKSVVCFConstants.SECOND_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, endCarrierSignal);
        refinedAttr.put(GATKSVVCFConstants.TOTAL_SPLIT_CARRIER_SIGNAL_ATTRIBUTE, totalCarrierSignal);
        refinedAttr.put(GATKSVVCFConstants.FIRST_SPLIT_POSITION_ATTRIBUTE, refinedFirstSite.getPosition());
        refinedAttr.put(GATKSVVCFConstants.SECOND_SPLIT_POSITION_ATTRIBUTE, refinedSecondSite.getPosition());

        final List<Genotype> genotypes = record.getGenotypes();
        final GenotypesContext newGenotypes = GenotypesContext.create(genotypes.size());
        for (final Genotype genotype : genotypes) {
            final String sample = genotype.getSampleName();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            genotypeBuilder.attribute(GATKSVVCFConstants.FIRST_SPLIT_READ_COUNT_ATTRIBUTE, refinedFirstSite.getCount(sample));
            genotypeBuilder.attribute(GATKSVVCFConstants.SECOND_SPLIT_READ_COUNT_ATTRIBUTE, refinedSecondSite.getCount(sample));
            newGenotypes.add(genotypeBuilder.make());
        }

        if (result.getPesrResult() != null) {
            final EvidenceStatUtils.PoissonTestResult discordantPairTest = result.getDiscordantPairTestResult().getTest();
            final Integer combinedCarrierSignal = EvidenceStatUtils.carrierSignalFraction(
                    discordantPairTest.getCarrierSignal() + bothsideResult.getCarrierSignal(),
                    discordantPairTest.getBackgroundSignal() + bothsideResult.getBackgroundSignal());
            refinedAttr.put(GATKSVVCFConstants.PESR_CARRIER_SIGNAL_ATTRIBUTE, combinedCarrierSignal);
            final Integer pesrQuality = Double.isNaN(result.getPesrResult().getP()) ?
                    null : EvidenceStatUtils.probToQual(result.getPesrResult().getP(), (byte) MAX_QUAL);
            refinedAttr.put(GATKSVVCFConstants.PESR_QUALITY_ATTRIBUTE, pesrQuality);
        }

        // Create new record
        return new SVCallRecord(record.getId(), record.getContigA(), newStart,
                refinedFirstSite.getStrand(), record.getContigB(), newEnd, refinedSecondSite.getStrand(),
                record.getType(), record.getComplexSubtype(), length, record.getAlgorithms(), record.getAlleles(),
                newGenotypes, refinedAttr, record.getFilters(), record.getLog10PError(), dictionary);
    }

    private static EvidenceStatUtils.PoissonTestResult calculateBothsideTest(final SplitReadSite startSite,
                                                                             final SplitReadSite endSite,
                                                                             final Set<String> carrierSamples,
                                                                             final Set<String> backgroundSamples,
                                                                             final Map<String, Double> sampleCoverageMap,
                                                                             final double representativeDepth) {
        final Map<String, Integer> sampleCountSums = new HashMap<>(SVUtils.hashMapCapacity(carrierSamples.size() + backgroundSamples.size()));
        for (final String sample : Sets.union(carrierSamples, backgroundSamples)) {
            sampleCountSums.put(sample, startSite.getCount(sample) + endSite.getCount(sample));
        }
        return EvidenceStatUtils.calculateOneSamplePoissonTest(sampleCountSums,
               carrierSamples, backgroundSamples, sampleCoverageMap, representativeDepth);
    }

    private static EvidenceStatUtils.PoissonTestResult calculatePESRTest(final SplitReadSite startSite,
                                                                         final SplitReadSite endSite,
                                                                         final DiscordantPairEvidenceTester.DiscordantPairTestResult discordantPairTestResult,
                                                                         final Set<String> carrierSamples,
                                                                         final Set<String> backgroundSamples,
                                                                         final Map<String, Double> sampleCoverageMap,
                                                                         final double representativeDepth) {
        final Map<String, Integer> sampleCountSums = new HashMap<>(SVUtils.hashMapCapacity(carrierSamples.size() + backgroundSamples.size()));
        final Map<String, Integer> discordantPairCounts = discordantPairTestResult.getSampleCounts();
        for (final String sample : Sets.union(carrierSamples, backgroundSamples)) {
            sampleCountSums.put(sample, startSite.getCount(sample) + endSite.getCount(sample) + discordantPairCounts.getOrDefault(sample, 0));
        }
        return EvidenceStatUtils.calculateOneSamplePoissonTest(sampleCountSums,
                carrierSamples, backgroundSamples, sampleCoverageMap, representativeDepth);
    }

    /**
     * Filters sites with position less than lower-bound
     *
     * @param evidence
     * @param lowerBound min position
     * @return filtered set of sites
     */
    private static List<SplitReadEvidence> filterSplitReadSitesLowerBound(final List<SplitReadEvidence> evidence, final int lowerBound) {
        return evidence.stream().filter(s -> s.getStart() >= lowerBound).collect(Collectors.toList());
    }

    /**
     * Filters sites with position greater than upper-bound
     *
     * @param evidence
     * @param upperBound min position
     * @return filtered set of sites
     */
    private static List<SplitReadEvidence> filterSplitReadSitesUpperBound(final List<SplitReadEvidence> evidence, final int upperBound) {
        return evidence.stream().filter(s -> s.getStart() <= upperBound).collect(Collectors.toList());
    }

    /**
     * Determines lower-bound on end site position (inclusive). For inter-chromosomal variants, boundaries are at the
     * start of the chromsome (any position is valid). For INS, {@link BreakpointRefiner#maxSplitReadCrossDistance}
     * is used to determine how far past the original breakpoint it can be. Otherwise, we just use the new start position.
     *
     * @param call
     * @param refinedStartPosition new start position of call
     * @return position
     */
    private int getEndLowerBound(final SVCallRecord call, final int refinedStartPosition) {
        if (!call.isIntrachromosomal()) {
            return 1;
        }
        if (call.getType().equals(GATKSVVCFConstants.StructuralVariantAnnotationType.INS)) {
            return refinedStartPosition - maxSplitReadCrossDistance;
        }
        return refinedStartPosition + 1;
    }

    /**
     * Same as {@link BreakpointRefiner#getEndLowerBound} but for upper-bound on start position.
     *
     * @param call
     * @param refinedEndPosition new end position of call
     * @return position
     */
    private int getStartUpperBound(final SVCallRecord call, final int refinedEndPosition) {
        if (!call.isIntrachromosomal()) {
            return Integer.MAX_VALUE;
        }
        if (call.getType().equals(GATKSVVCFConstants.StructuralVariantAnnotationType.INS)) {
            return refinedEndPosition + maxSplitReadCrossDistance;
        }
        return refinedEndPosition - 1;
    }

    public final class RefineResult {
        private final SplitReadSite first;
        private final SplitReadSite second;
        private final EvidenceStatUtils.PoissonTestResult bothsidesResult;
        private final DiscordantPairEvidenceTester.DiscordantPairTestResult discordantPairTestResult;
        private final EvidenceStatUtils.PoissonTestResult pesrResult;

        public RefineResult(final SplitReadSite first, final SplitReadSite second,
                            final EvidenceStatUtils.PoissonTestResult bothsidesResult,
                            final DiscordantPairEvidenceTester.DiscordantPairTestResult discordantPairTestResult,
                            final EvidenceStatUtils.PoissonTestResult pesrResult) {
            this.first = first;
            this.second = second;
            this.bothsidesResult = bothsidesResult;
            this.discordantPairTestResult = discordantPairTestResult;
            this.pesrResult = pesrResult;
        }

        public SplitReadSite getFirst() {
            return first;
        }

        public SplitReadSite getSecond() {
            return second;
        }

        public EvidenceStatUtils.PoissonTestResult getBothsidesResult() {
            return bothsidesResult;
        }

        public EvidenceStatUtils.PoissonTestResult getPesrResult() {
            return pesrResult;
        }

        public DiscordantPairEvidenceTester.DiscordantPairTestResult getDiscordantPairTestResult() {
            return discordantPairTestResult;
        }
    }
}
