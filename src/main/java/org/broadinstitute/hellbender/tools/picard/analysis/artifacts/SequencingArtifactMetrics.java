package org.broadinstitute.hellbender.tools.picard.analysis.artifacts;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.QualityUtil;

public final class SequencingArtifactMetrics {
    public static final String PRE_ADAPTER_SUMMARY_EXT = ".pre_adapter_summary_metrics";
    public static final String PRE_ADAPTER_DETAILS_EXT = ".pre_adapter_detail_metrics";
    public static final String BAIT_BIAS_SUMMARY_EXT = ".bait_bias_summary_metrics";
    public static final String BAIT_BIAS_DETAILS_EXT = ".bait_bias_detail_metrics";

    private static final double MIN_ERROR = 1e-10; // minimum error rate to report

    /**
     * Summary analysis of a single pre-adapter artifact.
     *
     * These artifacts occur on the original template strand, before the addition of adapters,
     * so they correlate with read number / orientation in a specific way.
     *
     * For example, the well-known "Oxo-G" artifact occurs when a G on the template
     * strand is oxidized, giving it an affinity for binding to A rather than the usual C.
     * Thus PCR will introduce apparent G>T substitutions in read 1 and C>A in read 2.
     * In the resulting alignments, a given G>T or C>A observation could either be:
     *
     * 1. a true mutation
     * 2. an OxoG artifact
     * 3. some other kind of artifact
     *
     * On average, we assume that 1 and 3 will not display this read number / orientation bias, so
     * their contributions will cancel out in the calculation.
     *
     */
    public static class PreAdapterSummaryMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;
        /** The (upper-case) original base on the reference strand. */
        public char REF_BASE;
        /** The (upper-case) alternative base that is called as a result of DNA damage. */
        public char ALT_BASE;

        /**
         * The total Phred-scaled Q-score for this artifact. A lower Q-score
         * means a higher probability that a REF_BASE:ALT_BASE observation
         * randomly picked from the data will be due to this artifact, rather
         * than a true variant.
         */
        public double TOTAL_QSCORE;

        /** The sequence context (reference bases surrounding the locus of interest) having the lowest Q-score among all contexts for this artifact. */
        public String WORST_CXT;
        /** The Q-score for the worst context. */
        public double WORST_CXT_QSCORE;

        /** The pre-context (reference bases leading up to the locus of interest) with the lowest Q-score. */
        public String WORST_PRE_CXT;
        /** The Q-score for the worst pre-context. */
        public double WORST_PRE_CXT_QSCORE;

        /** The post-context (reference bases trailing after the locus of interest) with the lowest Q-score. */
        public String WORST_POST_CXT;
        /** The Q-score for the worst post-context. */
        public double WORST_POST_CXT_QSCORE;

        /** A "nickname" of this artifact, if it is a known error mode. */
        public String ARTIFACT_NAME;

        /**
         * Label the artifacts corresponding to known error modes.
         */
        public void inferArtifactName() {
            if (this.REF_BASE == 'G' && this.ALT_BASE == 'T') this.ARTIFACT_NAME = "OxoG";
            else if (this.REF_BASE == 'C' && this.ALT_BASE == 'T') this.ARTIFACT_NAME = "Deamination";
            else this.ARTIFACT_NAME = "NA";
        }
    }

    /**
     * Summary analysis of a single bait bias artifact, also known as a reference bias artifact.
     *
     * These artifacts occur during or after the target selection step, and correlate with substitution
     * rates that are "biased", or higher for sites having one base on the reference/positive strand
     * relative to sites having the complementary base on that strand.
     *
     * For example, a G>T artifact during the target selection step might result in a higher
     * G>T / C>A substitution rate at sites with a G on the positive strand (and C on the negative),
     * relative to sites with the flip (C positive / G negative). This is known as the "G-Ref" artifact.
     *
     */
    public static class BaitBiasSummaryMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;
        /** The (upper-case) original base on the reference strand. */
        public char REF_BASE;
        /** The (upper-case) alternative base that is called as a result of DNA damage. */
        public char ALT_BASE;

        /**
         * The total Phred-scaled Q-score for this artifact. A lower Q-score
         * means a higher probability that a REF_BASE:ALT_BASE observation
         * randomly picked from the data will be due to this artifact, rather
         * than a true variant.
         */
        public double TOTAL_QSCORE;

        /** The sequence context (reference bases surrounding the locus of interest) having the lowest Q-score among all contexts for this artifact. */
        public String WORST_CXT;
        /** The Q-score for the worst context. */
        public double WORST_CXT_QSCORE;

        /** The pre-context (reference bases leading up to the locus of interest) with the lowest Q-score. */
        public String WORST_PRE_CXT;
        /** The Q-score for the worst pre-context. */
        public double WORST_PRE_CXT_QSCORE;

        /** The post-context (reference bases trailing after the locus of interest) with the lowest Q-score. */
        public String WORST_POST_CXT;
        /** The Q-score for the worst post-context. */
        public double WORST_POST_CXT_QSCORE;

        /** A "nickname" of this artifact, if it is a known error mode. */
        public String ARTIFACT_NAME;

        /**
         * Label the artifacts corresponding to known error modes.
         */
        public void inferArtifactName() {
            if (this.REF_BASE == 'G' && this.ALT_BASE == 'T') this.ARTIFACT_NAME = "Gref";
            else if (this.REF_BASE == 'C' && this.ALT_BASE == 'A') this.ARTIFACT_NAME = "Cref";
            else this.ARTIFACT_NAME = "NA";
        }
    }

    /**
     * Pre-adapter artifacts broken down by context.
     */
    public static class PreAdapterDetailMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;
        /** The (upper-case) original base on the reference strand. */
        public char REF_BASE;
        /** The (upper-case) alternative base that is called as a result of DNA damage. */
        public char ALT_BASE;

        /** The sequence context to which the analysis is constrained. */
        public String CONTEXT;

        /** The number of REF_BASE:REF_BASE alignments having a read number and orientation that supports the presence of this artifact. */
        public long PRO_REF_BASES;
        /** The number of REF_BASE:ALT_BASE alignments having a read number and orientation that supports the presence of this artifact. */
        public long PRO_ALT_BASES;
        /** The number of REF_BASE:REF_BASE alignments having a read number and orientation that refutes the presence of this artifact. */
        public long CON_REF_BASES;
        /** The number of REF_BASE:ALT_BASE alignments having a read number and orientation that refutes the presence of this artifact. */
        public long CON_ALT_BASES;

        /**
         * The estimated error rate due to this artifact.
         * Calculated as max(1e-10, (PRO_ALT_BASES - CON_ALT_BASES) / (PRO_ALT_BASES + PRO_REF_BASES + CON_ALT_BASES + CON_REF_BASES)).
         */
        public double ERROR_RATE;
        /** The Phred-scaled quality score of the artifact, calculated as -10 * log10(ERROR_RATE). */
        public double QSCORE;

        /**
         * Calculate the error rate given the raw counts. Negative rates are set to MIN_ERROR.
         */
        public void calculateDerivedStatistics() {
            this.ERROR_RATE = MIN_ERROR;
            final long totalBases = this.PRO_REF_BASES + this.PRO_ALT_BASES + this.CON_REF_BASES + this.CON_ALT_BASES;
            if (totalBases > 0) {
                final double rawErrorRate = (this.PRO_ALT_BASES - this.CON_ALT_BASES) / (double) totalBases;
                this.ERROR_RATE = Math.max(MIN_ERROR, rawErrorRate);
            }
            this.QSCORE = QualityUtil.getPhredScoreFromErrorProbability(this.ERROR_RATE);
        }
    }

    /**
     * Bait bias artifacts broken down by context.
     */
    public static class BaitBiasDetailMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;
        /** The (upper-case) original base on the reference strand. */
        public char REF_BASE;
        /** The (upper-case) alternative base that is called as a result of DNA damage. */
        public char ALT_BASE;

        /** The sequence context to which the analysis is constrained. */
        public String CONTEXT;

        /** The number of REF_BASE:REF_BASE alignments at sites with the given reference context. */
        public long FWD_CXT_REF_BASES;
        /** The number of REF_BASE:ALT_BASE alignments at sites with the given reference context. */
        public long FWD_CXT_ALT_BASES;
        /** The number of ~REF_BASE:~REF_BASE alignments at sites complementary to the given reference context. */
        public long REV_CXT_REF_BASES;
        /** The number of ~REF_BASE:~ALT_BASE alignments at sites complementary to the given reference context. */
        public long REV_CXT_ALT_BASES;

        /** The substitution rate of REF_BASE:ALT_BASE, calculated as max(1e-10, FWD_CXT_ALT_BASES / (FWD_CXT_ALT_BASES + FWD_CXT_REF_BASES)). */
        public double FWD_ERROR_RATE;
        /** The substitution rate of ~REF_BASE:~ALT_BASE, calculated as max(1e-10, REV_CXT_ALT_BASES / (REV_CXT_ALT_BASES + REV_CXT_REF_BASES)). */
        public double REV_ERROR_RATE;

        /**
         * The bait bias error rate, calculated as max(1e-10, FWD_ERROR_RATE - REV_ERROR_RATE).
         */
        public double ERROR_RATE;
        /** The Phred-scaled quality score of the artifact, calculated as -10 * log10(ERROR_RATE). */
        public double QSCORE;

        /**
         * Calculate the error rate given the raw counts. Negative rates are set to MIN_ERROR.
         */
        public void calculateDerivedStatistics() {
            this.FWD_ERROR_RATE = MIN_ERROR;
            final long fwdBases = this.FWD_CXT_REF_BASES + this.FWD_CXT_ALT_BASES;
            if (fwdBases > 0) {
                final double fwdErr = this.FWD_CXT_ALT_BASES / (double) fwdBases;
                this.FWD_ERROR_RATE = Math.max(MIN_ERROR, fwdErr);
            }

            this.REV_ERROR_RATE = MIN_ERROR;
            final long revBases = this.REV_CXT_REF_BASES + this.REV_CXT_ALT_BASES;
            if (revBases > 0) {
                final double revErr = this.REV_CXT_ALT_BASES / (double) revBases;
                this.REV_ERROR_RATE = Math.max(MIN_ERROR, revErr);
            }

            this.ERROR_RATE = Math.max(MIN_ERROR, this.FWD_ERROR_RATE - this.REV_ERROR_RATE);
            this.QSCORE = QualityUtil.getPhredScoreFromErrorProbability(this.ERROR_RATE);
        }
    }

    /**
     * Little container for passing around a detail metrics pair.
     */
    static class DetailPair {
        final PreAdapterDetailMetrics preAdapterMetrics;
        final BaitBiasDetailMetrics baitBiasMetrics;

        DetailPair(final PreAdapterDetailMetrics preAdapterMetrics, final BaitBiasDetailMetrics baitBiasMetrics) {
            this.preAdapterMetrics = preAdapterMetrics;
            this.baitBiasMetrics = baitBiasMetrics;
        }
    }

    /**
     * Little container for passing around a summary metrics pair.
     */
    static class SummaryPair {
        final PreAdapterSummaryMetrics preAdapterMetrics;
        final BaitBiasSummaryMetrics baitBiasMetrics;

        SummaryPair(final PreAdapterSummaryMetrics preAdapterMetrics, final BaitBiasSummaryMetrics baitBiasMetrics) {
            this.preAdapterMetrics = preAdapterMetrics;
            this.baitBiasMetrics = baitBiasMetrics;
        }
    }
}
