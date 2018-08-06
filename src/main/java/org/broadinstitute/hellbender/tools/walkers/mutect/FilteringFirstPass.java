package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class FilteringFirstPass {
    final List<FilterResult> filterResults;
    final Map<String, ImmutablePair<String, Integer>> filteredPhasedCalls;

    public FilteringFirstPass() {
        filterResults = new ArrayList<>();
        filteredPhasedCalls = new HashMap<>();
    }

    public void add(final FilterResult filterResult, final VariantContext vc) {
        filterResults.add(filterResult);

        if (!filterResult.getFilters().isEmpty() && hasPhaseInfo(vc)) {
            final String pgt = vc.getAttributeAsString(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "");
            final String pid = vc.getAttributeAsString(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, "");
            final int position = vc.getStart();
            filteredPhasedCalls.put(pgt, new ImmutablePair<>(pid, position));
        }
    }

    public Mutect2FilterSummary calculateFilterStats(final double requestedFPR){
        final Mutect2FilterSummary filterSummary = new Mutect2FilterSummary();

        final double[] readOrientationPosteriors = getFilterResults().stream()
                .filter(r -> r.getFilters().isEmpty())
                .mapToDouble(r -> r.getReadOrientationPosterior())
                .toArray();

        final Mutect2FilterSummary.FilterStats readOrientationFilterStats = calculateThresholdForReadOrientationFilter(readOrientationPosteriors, requestedFPR);
        filterSummary.addNewFilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, readOrientationFilterStats);

        return filterSummary;
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
    public static Mutect2FilterSummary.FilterStats calculateThresholdForReadOrientationFilter(final double[] posteriors, final double requestedFPR){
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
                        new Mutect2FilterSummary.FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, posteriors[i-1],
                                cumulativeExpectedFPs, i-1, cumulativeExpectedFPs/i, requestedFPR) :
                        new Mutect2FilterSummary.FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, thresholdForFilteringAll,
                                0.0, 0, 0.0, requestedFPR);
            }

            cumulativeExpectedFPs += posterior;
        }

        // If the expected FP rate never exceeded the max tolerable value, then we can let everything pass
        return new Mutect2FilterSummary.FilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME, thresholdForFilteringNone,
                cumulativeExpectedFPs, numPassingVariants, cumulativeExpectedFPs/numPassingVariants, requestedFPR);
    }

    private static boolean hasPhaseInfo(VariantContext vc) {
        return vc.hasAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) && vc.hasAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY);
    }

    public List<FilterResult> getFilterResults() {
        return filterResults;
    }



}
