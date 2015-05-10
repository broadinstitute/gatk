package org.broadinstitute.hellbender.tools.picard.analysis.artifacts;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ListMap;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.SequencingArtifactMetrics.*;

import java.util.*;

/**
 * Keeps track of artifact counts, and extracts metrics once accumulation is finished.
 */
final class ArtifactCounter {
    private final String sampleAlias;
    private final String library;

    private final Set<String> fullContexts;
    private final Map<String, String> leadingContextMap;
    private final Map<String, String> trailingContextMap;
    private final Map<String, String> zeroContextMap;

    private final ContextAccumulator fullContextAccumulator;
    private final ContextAccumulator halfContextAccumulator;
    private final ContextAccumulator zeroContextAccumulator;

    private final List<PreAdapterSummaryMetrics> preAdapterSummaryMetricsList;
    private final List<PreAdapterDetailMetrics> preAdapterDetailMetricsList;
    private final List<BaitBiasSummaryMetrics> baitBiasSummaryMetricsList;
    private final List<BaitBiasDetailMetrics> baitBiasDetailMetricsList;

    public ArtifactCounter(final String sampleAlias, final String library, final int contextSize, final boolean expectedTandemReads) {
        this.sampleAlias = sampleAlias;
        this.library = library;

        // define the contexts
        this.fullContexts = new HashSet<>();
        for (final byte[] kmer : SequenceUtil.generateAllKmers(2 * contextSize + 1)) {
            this.fullContexts.add(StringUtil.bytesToString(kmer));
        }

        // the half contexts specify either leading or trailing bases. the zero context is just the center.
        // NB: we use N to represent a wildcard base, rather than an ambiguous base. It's assumed that all of the input
        // contexts are unambiguous, and that any actual N's in the data have been dealt with elsewhere.
        final String padding = StringUtil.repeatCharNTimes('N', contextSize);
        this.leadingContextMap = new HashMap<>();
        this.trailingContextMap = new HashMap<>();
        this.zeroContextMap = new HashMap<>();
        for (final String context : this.fullContexts) {
            final String leading = context.substring(0, contextSize);
            final String trailing = context.substring(contextSize + 1, context.length());
            final char center = context.charAt(contextSize);
            this.leadingContextMap.put(context, leading + center + padding);
            this.trailingContextMap.put(context, padding + center + trailing);
            this.zeroContextMap.put(context, padding + center + padding);
        }

        // set up the accumulators
        final Set<String> halfContexts = new HashSet<>();
        halfContexts.addAll(leadingContextMap.values());
        halfContexts.addAll(trailingContextMap.values());
        final Set<String> zeroContexts = new HashSet<>();
        zeroContexts.addAll(zeroContextMap.values());

        this.fullContextAccumulator = new ContextAccumulator(fullContexts, expectedTandemReads);
        this.halfContextAccumulator = new ContextAccumulator(halfContexts, expectedTandemReads);
        this.zeroContextAccumulator = new ContextAccumulator(zeroContexts, expectedTandemReads);

        // these will get populated in the final step
        preAdapterSummaryMetricsList = new ArrayList<>();
        preAdapterDetailMetricsList = new ArrayList<>();
        baitBiasSummaryMetricsList = new ArrayList<>();
        baitBiasDetailMetricsList = new ArrayList<>();
    }

    /**
     * Add a record to all the accumulators.
     */
    public void countRecord(final String refContext, final char calledBase, final SAMRecord rec) {
        this.fullContextAccumulator.countRecord(refContext, calledBase, rec);
        this.halfContextAccumulator.countRecord(this.leadingContextMap.get(refContext), calledBase, rec);
        this.halfContextAccumulator.countRecord(this.trailingContextMap.get(refContext), calledBase, rec);
        this.zeroContextAccumulator.countRecord(this.zeroContextMap.get(refContext), calledBase, rec);
    }

    /**
     * Stop counting, tally things up, and extract metrics.
     */
    public void finish() {
        final ListMap<Transition, DetailPair> allDetailMetrics = getDetailMetrics();
        final Map<Transition, SummaryPair> allSummaryMetrics = getSummaryMetrics();

        for (final Transition transition : Transition.altValues()) {
            final SummaryPair summary = allSummaryMetrics.get(transition);
            final List<DetailPair> details = allDetailMetrics.get(transition);
            preAdapterSummaryMetricsList.add(summary.preAdapterMetrics);
            baitBiasSummaryMetricsList.add(summary.baitBiasMetrics);
            for (final DetailPair detail : details) {
                preAdapterDetailMetricsList.add(detail.preAdapterMetrics);
                baitBiasDetailMetricsList.add(detail.baitBiasMetrics);
            }
        }
    }

    public List<PreAdapterSummaryMetrics> getPreAdapterSummaryMetrics() { return preAdapterSummaryMetricsList; }
    public List<PreAdapterDetailMetrics> getPreAdapterDetailMetrics() { return preAdapterDetailMetricsList; }
    public List<BaitBiasSummaryMetrics> getBaitBiasSummaryMetrics() { return baitBiasSummaryMetricsList; }
    public List<BaitBiasDetailMetrics> getBaitBiasDetailMetrics() { return baitBiasDetailMetricsList; }

    /**
     * Core method to compute summary metrics. For each transition, we report:
     * 1. the total Q-score across all contexts
     * 2. the worst full context and its Q-score
     * 3. the worst leading context and its Q-score
     * 4. the worst trailing context and its Q-score
     *
     */
    private Map<Transition, SummaryPair> getSummaryMetrics() {
        final Map<Transition, SummaryPair> summaryMetricsMap = new HashMap<>();

        // extract the detail metrics from each accumulator
        final ListMap<Transition, DetailPair> fullMetrics = this.fullContextAccumulator.calculateMetrics(sampleAlias, library);
        final ListMap<Transition, DetailPair> halfMetrics = this.halfContextAccumulator.calculateMetrics(sampleAlias, library);
        final ListMap<Transition, DetailPair> zeroMetrics = this.zeroContextAccumulator.calculateMetrics(sampleAlias, library);

        // compute the summary metrics - one row for each transition
        for (final Transition transition : Transition.altValues()) {
            final List<DetailPair> fullMetricsForTransition = fullMetrics.get(transition);
            final List<DetailPair> zeroMetricsForTransition = zeroMetrics.get(transition);
            if (zeroMetricsForTransition.size() != 1) {
                throw new GATKException("Should have exactly one context-free metric pair for transition: " + transition);
            }

            // we want to report on leading / trailing contexts separately
            final List<DetailPair> leadingMetricsForTransition = new ArrayList<>();
            final List<DetailPair> trailingMetricsForTransition = new ArrayList<>();
            for (final DetailPair metrics : halfMetrics.get(transition)) {
                // first make sure they're the same context
                if (!metrics.preAdapterMetrics.CONTEXT.equals(metrics.baitBiasMetrics.CONTEXT)) {
                    throw new GATKException("Input detail metrics are not matched up properly - contexts differ.");
                }
                final boolean isLeading = this.leadingContextMap.containsValue(metrics.preAdapterMetrics.CONTEXT);
                final boolean isTrailing = this.trailingContextMap.containsValue(metrics.preAdapterMetrics.CONTEXT);
                // if the original contextSize is 0, there's no difference between leading and trailing, so add it to both
                if (isLeading) leadingMetricsForTransition.add(metrics);
                if (isTrailing) trailingMetricsForTransition.add(metrics);
            }

            // get the worst cases
            final DetailPair totalMetric = zeroMetricsForTransition.get(0);
            final DetailPair worstFullMetric = getWorstMetrics(fullMetricsForTransition);
            final DetailPair worstLeadingMetric = getWorstMetrics(leadingMetricsForTransition);
            final DetailPair worstTrailingMetric = getWorstMetrics(trailingMetricsForTransition);

            // construct the actual summary metrics - a combination of all the data we've just extracted
            final PreAdapterSummaryMetrics preAdapterSummaryMetrics = new PreAdapterSummaryMetrics();
            final BaitBiasSummaryMetrics baitBiasSummaryMetrics = new BaitBiasSummaryMetrics();

            preAdapterSummaryMetrics.SAMPLE_ALIAS = this.sampleAlias;
            preAdapterSummaryMetrics.LIBRARY = this.library;
            preAdapterSummaryMetrics.REF_BASE = transition.ref();
            preAdapterSummaryMetrics.ALT_BASE = transition.call();
            preAdapterSummaryMetrics.TOTAL_QSCORE = totalMetric.preAdapterMetrics.QSCORE;
            preAdapterSummaryMetrics.WORST_CXT = worstFullMetric.preAdapterMetrics.CONTEXT;
            preAdapterSummaryMetrics.WORST_CXT_QSCORE = worstFullMetric.preAdapterMetrics.QSCORE;
            preAdapterSummaryMetrics.WORST_PRE_CXT = worstLeadingMetric.preAdapterMetrics.CONTEXT;
            preAdapterSummaryMetrics.WORST_PRE_CXT_QSCORE = worstLeadingMetric.preAdapterMetrics.QSCORE;
            preAdapterSummaryMetrics.WORST_POST_CXT = worstTrailingMetric.preAdapterMetrics.CONTEXT;
            preAdapterSummaryMetrics.WORST_POST_CXT_QSCORE = worstTrailingMetric.preAdapterMetrics.QSCORE;
            preAdapterSummaryMetrics.inferArtifactName();

            baitBiasSummaryMetrics.SAMPLE_ALIAS = this.sampleAlias;
            baitBiasSummaryMetrics.LIBRARY = this.library;
            baitBiasSummaryMetrics.REF_BASE = transition.ref();
            baitBiasSummaryMetrics.ALT_BASE = transition.call();
            baitBiasSummaryMetrics.TOTAL_QSCORE = totalMetric.baitBiasMetrics.QSCORE;
            baitBiasSummaryMetrics.WORST_CXT = worstFullMetric.baitBiasMetrics.CONTEXT;
            baitBiasSummaryMetrics.WORST_CXT_QSCORE = worstFullMetric.baitBiasMetrics.QSCORE;
            baitBiasSummaryMetrics.WORST_PRE_CXT = worstLeadingMetric.baitBiasMetrics.CONTEXT;
            baitBiasSummaryMetrics.WORST_PRE_CXT_QSCORE = worstLeadingMetric.baitBiasMetrics.QSCORE;
            baitBiasSummaryMetrics.WORST_POST_CXT = worstTrailingMetric.baitBiasMetrics.CONTEXT;
            baitBiasSummaryMetrics.WORST_POST_CXT_QSCORE = worstTrailingMetric.baitBiasMetrics.QSCORE;
            baitBiasSummaryMetrics.inferArtifactName();

            // add the finalized metrics to the map
            summaryMetricsMap.put(transition, new SummaryPair(preAdapterSummaryMetrics, baitBiasSummaryMetrics));
        }
        return summaryMetricsMap;
    }

    private ListMap<Transition, DetailPair> getDetailMetrics() {
        return this.fullContextAccumulator.calculateMetrics(this.sampleAlias, this.library);
    }

    /**
     * Given a list of detail metrics, get the worst pre-adapter metrics, and independently from that get the worst bait bias metrics
     * (in terms of Q-score).
     */
    private DetailPair getWorstMetrics(final List<DetailPair> metrics) {
        PreAdapterDetailMetrics worstPreAdapterMetrics = null;
        BaitBiasDetailMetrics worstBaitBiasMetrics = null;
        for (final DetailPair m : metrics) {
            if (worstPreAdapterMetrics == null || m.preAdapterMetrics.QSCORE < worstPreAdapterMetrics.QSCORE) worstPreAdapterMetrics = m.preAdapterMetrics;
            if (worstBaitBiasMetrics == null || m.baitBiasMetrics.QSCORE < worstBaitBiasMetrics.QSCORE) worstBaitBiasMetrics = m.baitBiasMetrics;
        }
        return new DetailPair(worstPreAdapterMetrics, worstBaitBiasMetrics);
    }
}
