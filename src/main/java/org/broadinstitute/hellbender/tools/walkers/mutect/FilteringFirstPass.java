package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Stores the results of the first pass of {@link FilterMutectCalls}, a purely online step in which each variant is
 * not "aware" of other variants, and learns various global properties necessary for a more refined second step.
 */
public class FilteringFirstPass {
    private final List<FilterResult> filterResults;

    // for each PID, the positions with PGTs of filtered genotypes
    private final Map<String, ImmutablePair<Integer, Set<String>>> filteredPhasedCalls;
    private final Map<String, FilterStats> filterStats;
    private final Set<String> normalSamples;
    boolean readyForSecondPass;

    public FilteringFirstPass(final Set<String> normalSamples) {
        filterResults = new ArrayList<>();
        filteredPhasedCalls = new HashMap<>();
        filterStats = new HashMap<>();
        readyForSecondPass = false;
        this.normalSamples = normalSamples;
    }

    public boolean isReadyForSecondPass() { return readyForSecondPass; }

    public FilterStats getFilterStats(final String filterName){
        Utils.validateArg(filterStats.containsKey(filterName), "invalid filter name: " + filterName);
        return filterStats.get(filterName);
    }

    public boolean isOnFilteredHaplotype(final VariantContext vc, final int maxDistance) {

        final Genotype tumorGenotype = getTumorGenotypeWithGreatestAlleleFraction(vc);

        if (!hasPhaseInfo(tumorGenotype)) {
            return false;
        }

        final String pgt = (String) tumorGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "");
        final String pid = (String) tumorGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, "");
        final int position = vc.getStart();

        final Pair<Integer, Set<String>> filteredCall = filteredPhasedCalls.get(pid);
        if (filteredCall == null) {
            return false;
        }

        // Check that vc occurs on a filtered haplotype
        return filteredCall.getRight().contains(pgt) && Math.abs(filteredCall.getLeft() - position) <= maxDistance;
    }

    private Genotype getTumorGenotypeWithGreatestAlleleFraction(final VariantContext vc) {
        return vc.getGenotypes().stream()
                    .filter(g ->  !normalSamples.contains(g.getSampleName()))
                    .max(Comparator.comparingDouble(g -> MathUtils.arrayMax(GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                            () -> new double[] {0.0}, 0.0)))).get();
    }

    public void add(final FilterResult filterResult, final VariantContext vc) {
        filterResults.add(filterResult);

        final int position = vc.getStart();
        if (!filterResult.getFilters().isEmpty()) {
            final Map<String, Set<String>> phasedGTsForEachPhaseID = vc.getGenotypes().stream()
                    .filter(gt -> !normalSamples.contains(gt.getSampleName()))
                    .filter(FilteringFirstPass::hasPhaseInfo)
                    .collect(Collectors.groupingBy(g -> (String) g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, ""),
                            Collectors.mapping(g -> (String) g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, ""), Collectors.toSet())));

            for (final Map.Entry<String, Set<String>> pidAndPgts : phasedGTsForEachPhaseID.entrySet()) {
                filteredPhasedCalls.put(pidAndPgts.getKey(), new ImmutablePair<>(position, pidAndPgts.getValue()));
            }
        }
    }

    public void learnModelForSecondPass(final double requestedFPR){
        final double[] readOrientationPosteriors = getFilterResults().stream()
                .filter(r -> r.getFilters().isEmpty())
                .mapToDouble(r -> r.getReadOrientationPosterior())
                .toArray();

        final FilterStats readOrientationFilterStats = calculateThresholdForReadOrientationFilter(readOrientationPosteriors, requestedFPR);
        filterStats.put(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, readOrientationFilterStats);
        readyForSecondPass = true;
    }

    /**
     *
     * Compute the filtering threshold that ensures that the false positive rate among the resulting pass variants
     * will not exceed the requested false positive rate
     *
     * @param posteriors A list of posterior probabilities, which gets sorted
     * @param requestedFPR We set the filtering threshold such that the FPR doesn't exceed this value
     * @return
     */
    public static FilterStats calculateThresholdForReadOrientationFilter(final double[] posteriors, final double requestedFPR){
        ParamUtils.isPositiveOrZero(requestedFPR, "requested FPR must be non-negative");
        final double thresholdForFilteringNone = 1.0;
        final double thresholdForFilteringAll = 0.0;

        Arrays.sort(posteriors);

        final int numPassingVariants = posteriors.length;
        double cumulativeExpectedFPs = 0.0;

        for (int i = 0; i < numPassingVariants; i++){
            final double posterior = posteriors[i];

            // One can show that the cumulative error rate is monotonically increasing in i
            final double expectedFPR = (cumulativeExpectedFPs + posterior) / (i + 1);
            if (expectedFPR > requestedFPR){
                return i > 0 ?
                        new FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, posteriors[i-1],
                                cumulativeExpectedFPs, i-1, cumulativeExpectedFPs/i, requestedFPR) :
                        new FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, thresholdForFilteringAll,
                                0.0, 0, 0.0, requestedFPR);
            }

            cumulativeExpectedFPs += posterior;
        }

        // If the expected FP rate never exceeded the max tolerable value, then we can let everything pass
        return new FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, thresholdForFilteringNone,
                cumulativeExpectedFPs, numPassingVariants, cumulativeExpectedFPs/numPassingVariants, requestedFPR);
    }

    public static boolean hasPhaseInfo(final Genotype genotype) {
        return genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) && genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY);
    }

    public List<FilterResult> getFilterResults() {
        return filterResults;
    }

    public static class FilterStats {
        private final String filterName;
        private final double threshold;
        private final double expectedNumFPs;
        private final int numPassingVariants;
        private final double expectedFPR;
        private final double requestedFPR;

        public FilterStats(final String filterName, final double threshold, final double expectedNumFPs,
                           final int numPassingVariants, final double expectedFPR, final double requestedFPR){
            this.filterName = filterName;
            this.threshold = threshold;
            this.expectedNumFPs = expectedNumFPs;
            this.numPassingVariants = numPassingVariants;
            this.expectedFPR = expectedFPR;
            this.requestedFPR = requestedFPR;
        }

        public String getFilterName() { return filterName; }

        public double getExpectedNumFPs() { return expectedNumFPs; }

        public int getNumPassingVariants() { return numPassingVariants; }

        public double getThreshold() { return threshold; }

        public double getExpectedFPR() { return expectedFPR; }

        public double getRequestedFPR() { return requestedFPR; }

    }

    private enum M2FilterStatsTableColumn {
        FILTER_NAME("filter_name"),
        THRESHOLD("threshold"),
        EXPECTED_FALSE_POSITIVES("expected_fps"),
        EXPECTED_FALSE_POSITIVE_RATE("expected_fpr"),
        REQUESTED_FALSE_POSITIVE_RATE("requested_fpr"),
        NUM_PASSING_VARIANTS("num_passing_variants");

        private String columnName;

        M2FilterStatsTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() { return columnName; }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static class Mutect2FilterStatsWriter extends TableWriter<FilterStats> {
        private Mutect2FilterStatsWriter(final File output) throws IOException {
            super(output, M2FilterStatsTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final FilterStats stats, final DataLine dataLine) {
            dataLine.set(M2FilterStatsTableColumn.FILTER_NAME.toString(), stats.getFilterName())
                    .set(M2FilterStatsTableColumn.THRESHOLD.toString(), stats.getThreshold())
                    .set(M2FilterStatsTableColumn.EXPECTED_FALSE_POSITIVES.toString(), stats.getExpectedNumFPs())
                    .set(M2FilterStatsTableColumn.EXPECTED_FALSE_POSITIVE_RATE.toString(), stats.getExpectedFPR())
                    .set(M2FilterStatsTableColumn.REQUESTED_FALSE_POSITIVE_RATE.toString(), stats.getRequestedFPR())
                    .set(M2FilterStatsTableColumn.NUM_PASSING_VARIANTS.toString(), stats.getNumPassingVariants());
        }
    }

    public void writeM2FilterSummary(final File outputTable) {
        try (Mutect2FilterStatsWriter writer = new Mutect2FilterStatsWriter(outputTable)) {
            writer.writeAllRecords(filterStats.values());
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable), e);
        }
    }
}
