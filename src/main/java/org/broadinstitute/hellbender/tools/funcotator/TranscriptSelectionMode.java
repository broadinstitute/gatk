package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;

import java.util.Comparator;
import java.util.Set;

/**
 * The manner to select a single transcript from a set of transcripts to report as the "best" or main transcript.
 * That is, the transcript that has the most information on it in the output.
 */
public enum TranscriptSelectionMode {

    /**
     * BEST_EFFECT
     *
     * Select a transcript to be reported with details with priority on effect according to the folowing list of selection criteria:
     *
     * Choose the transcript that is on the custom list specified by the user. If no list was specified, treat as if no transcripts were on the list (tie).
     * In case of tie, choose the transcript that yields the variant classification highest on the variant classification rank list (see below).
     * If still a tie, choose the transcript with highest level of curation. Note that this means lower number is better for level (see below).
     * If still a tie, choose the transcript with the best appris annotation (see below).
     * If still a tie, choose the transcript with the longest transcript sequence length.
     * If still a tie, choose the first transcript, alphabetically.
     *
     * Levels of Curation:
     *     1 (verified loci)
     *     2 manually annotated loci
     *     3 automatically annotated loci
     *
     * Variant Classification Scores (See {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} as well):
     *      De_novo_Start_OutOfFrame    0
     *      Nonsense_Mutation           0
     *      Nonstop_Mutation            0
     *      Missense_Mutation           1
     *      De_novo_Start_InFrame       1
     *      In_Frame_Del                1
     *      In_Frame_Ins                1
     *      Frame_Shift_Del             2
     *      Frame_Shift_Ins             2
     *      Frame_Shift_Sub             2
     *      Start_Codon_SNP             3
     *      Start_Codon_Del             3
     *      Start_Codon_Ins             3
     *      Start_Codon_DNP             3
     *      Start_Codon_TNP             3
     *      Start_Codon_ONP             3
     *      Stop_Codon_SNP              3
     *      Stop_Codon_Del              3
     *      Stop_Codon_Ins              3
     *      Stop_Codon_DNP              3
     *      Stop_Codon_TNP              3
     *      Stop_Codon_ONP              3
     *      Splice_Site                 4
     *      Splice_Site_SNP             4
     *      Splice_Site_Del             4
     *      Splice_Site_Ins             4
     *      Splice_Site_DNP             4
     *      Splice_Site_TNP             4
     *      Splice_Site_ONP             4
     *      Splice_Site                 4
     *      miRNA                       4
     *      Silent                      5
     *      3UTR                        6
     *      5UTR                        6
     *      Intron                      7
     *      5Flank                      8
     *      3Flank                      8
     *      Non-coding_Transcript       9
     *      IGR	                        20
     *      TX-REF-MISMATCH             100
     *
     * APPRIS Ranks (http://appris.bioinfo.cnio.es/):
     *
     *      appris_principal
     *      appris_candidate_highest_score
     *      appris_candidate_longest_ccds
     *      appris_candidate_ccds
     *      appris_candidate_longest_seq
     *      appris_candidate_longest
     *      appris_candidate
     *      no appris tag present
     */
    BEST_EFFECT {
        public Comparator<GencodeFuncotation> getComparator(final Set<String> userRequestedTranscripts) {
            return new BestEffectGencodeFuncotationComparator(userRequestedTranscripts);
        }
    },

    /**
     * CANONICAL
     *
     * Select a transcript to be reported with details with priority on canonical order according to the folowing list of selection criteria:
     *
     * Choose the transcript that is on the custom list specified by the user. If no list was specified, treat as if all transcripts were on the list (tie).
     * In case of tie, choose the transcript with highest level of curation. Note that this means lower number is better for level (see below).
     * If still a tie, choose the transcript that yields the variant classification highest on the variant classification rank list (see below).
     * If still a tie, choose the transcript with the best appris annotation (see below).
     * If still a tie, choose the transcript with the longest transcript sequence length.
     * If still a tie, choose the first transcript, alphabetically.
     *
     * Levels of Curation:
     *     1 (verified loci)
     *     2 manually annotated loci
     *     3 automatically annotated loci
     *
     * Variant Classification Scores (See {@link GencodeFuncotation.VariantClassification} as well):
     *      De_novo_Start_OutOfFrame    0
     *      Nonsense_Mutation           0
     *      Nonstop_Mutation            0
     *      Missense_Mutation           1
     *      De_novo_Start_InFrame       1
     *      In_Frame_Del                1
     *      In_Frame_Ins                1
     *      Frame_Shift_Del             2
     *      Frame_Shift_Ins             2
     *      Frame_Shift_Sub             2
     *      Start_Codon_SNP             3
     *      Start_Codon_Del             3
     *      Start_Codon_Ins             3
     *      Start_Codon_DNP             3
     *      Start_Codon_TNP             3
     *      Start_Codon_ONP             3
     *      Stop_Codon_SNP              3
     *      Stop_Codon_Del              3
     *      Stop_Codon_Ins              3
     *      Stop_Codon_DNP              3
     *      Stop_Codon_TNP              3
     *      Stop_Codon_ONP              3
     *      Splice_Site                 4
     *      Splice_Site_SNP             4
     *      Splice_Site_Del             4
     *      Splice_Site_Ins             4
     *      Splice_Site_DNP             4
     *      Splice_Site_TNP             4
     *      Splice_Site_ONP             4
     *      Splice_Site                 4
     *      miRNA                       4
     *      Silent                      5
     *      3UTR                        6
     *      5UTR                        6
     *      Intron                      7
     *      5Flank                      8
     *      3Flank                      8
     *      Non-coding_Transcript       9
     *      IGR	                        20
     *      TX-REF-MISMATCH             100
     *
     * APPRIS Ranks (http://appris.bioinfo.cnio.es/):
     *
     *      appris_principal
     *      appris_candidate_highest_score
     *      appris_candidate_longest_ccds
     *      appris_candidate_ccds
     *      appris_candidate_longest_seq
     *      appris_candidate_longest
     *      appris_candidate
     *      no appris tag present
     */
    CANONICAL {
        public Comparator<GencodeFuncotation> getComparator(final Set<String> userRequestedTranscripts) {
            return new CanonicalGencodeFuncotationComparator(userRequestedTranscripts);
        }
    },

    /**
     * Same as CANONICAL, but indicates that no transcripts should be dropped.  Render all overlapping transcripts.
     */
    ALL {
        public Comparator<GencodeFuncotation> getComparator(final Set<String> userRequestedTranscripts) {
            return new CanonicalGencodeFuncotationComparator(userRequestedTranscripts);
        }
    };

    public abstract Comparator<GencodeFuncotation> getComparator(final Set<String> userRequestedTranscripts);

    private static class ComparatorByUserTranscript implements Comparator<GencodeFuncotation> {

        private final Set<String> userRequestedTranscripts;

        public ComparatorByUserTranscript( final Set<String> userRequestedTranscripts ) {
            this.userRequestedTranscripts = userRequestedTranscripts;
        }

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // Choose the transcript that is on the custom list specified by the user:
            if ( FuncotatorUtils.isFuncotationInTranscriptList(a, userRequestedTranscripts) && (!FuncotatorUtils.isFuncotationInTranscriptList(b, userRequestedTranscripts)) ) {
                return -1;
            }
            else if ( (!FuncotatorUtils.isFuncotationInTranscriptList(a, userRequestedTranscripts)) && FuncotatorUtils.isFuncotationInTranscriptList(b, userRequestedTranscripts) ) {
                return 1;
            }
            else {
                return 0;
            }
        }
    }

    private static class ComparatorByIgrStatus implements Comparator<GencodeFuncotation> {

        public ComparatorByIgrStatus(){}

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // Check to see if one is an IGR.  IGR's have only a subset of the information in them, so it's easier to
            // order them if they're IGRs:
            if ( (b.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) &&
                    (!a.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) ) {
                return -1;
            }
            else if ( (a.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) &&
                    (!b.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) ) {
                return 1;
            }
            else {
                return 0;
            }
        }
    }

    private static class ComparatorByVariantClassification implements Comparator<GencodeFuncotation> {

        public ComparatorByVariantClassification(){}

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // Check highest variant classification:
            if ( a.getVariantClassification().getSeverity() < b.getVariantClassification().getSeverity() ) {
                return -1;
            }
            else if ( a.getVariantClassification().getSeverity() > b.getVariantClassification().getSeverity() ) {
                return 1;
            }
            else {
                return 0;
            }
        }
    }

    private static class ComparatorByProteinCodingStatus implements Comparator<GencodeFuncotation> {

        public ComparatorByProteinCodingStatus(){}

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // Is it protein coding?
            final boolean isAProteinCoding = a.getGeneTranscriptType() == GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING;
            final boolean isBProteinCoding = b.getGeneTranscriptType() == GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING;
            if ( isAProteinCoding != isBProteinCoding ) {
                if ( isAProteinCoding ) {
                    return -1;
                }
                else {
                    return 1;
                }
            }
            else {
                return 0;
            }
        }
    }

    private static class ComparatorByLocusLevel implements Comparator<GencodeFuncotation> {

        public ComparatorByLocusLevel(){}

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // Check locus/curation levels:
            if ( (a.getLocusLevel() != null) && (b.getLocusLevel() == null) ) {
                return -1;
            }
            else if ( (a.getLocusLevel() == null ) && (b.getLocusLevel() != null) ) {
                return 1;
            }
            else if ( (a.getLocusLevel() == null ) && (b.getLocusLevel() == null) ) {
                return 0;
            }
            else {
                return a.getLocusLevel().compareTo( b.getLocusLevel() );
            }
        }
    }

    private static class ComparatorByApprisRank implements Comparator<GencodeFuncotation> {

        public ComparatorByApprisRank(){}

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // Check the appris annotation:
            if ( (a.getApprisRank() != null) && (b.getApprisRank() == null) ) {
                return -1;
            }
            else if ( (a.getApprisRank() == null ) && (b.getApprisRank() != null) ) {
                return 1;
            }
            else if ( (a.getApprisRank() != null) && (b.getApprisRank() != null) && (!a.getApprisRank().equals(b.getApprisRank())) ) {
                return a.getApprisRank().compareTo( b.getApprisRank() );
            }
            else {
                return 0;
            }
        }
    }

    private static class ComparatorByTranscriptSequenceLength implements Comparator<GencodeFuncotation> {

        public ComparatorByTranscriptSequenceLength(){}

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // Check transcript sequence length:
            // Note that since we want longer transcripts to be sorted earlier in the list, we reverse b and a in the last check.
            if ( (a.getTranscriptLength() != null) && (b.getTranscriptLength() == null) ) {
                return -1;
            }
            else if ( (a.getTranscriptLength() == null ) && (b.getTranscriptLength() != null) ) {
                return 1;
            }
            else if ( (a.getTranscriptLength() != null) && (b.getTranscriptLength() != null) && (!a.getTranscriptLength().equals(b.getTranscriptLength())) ) {
                return b.getTranscriptLength().compareTo( a.getTranscriptLength() );
            }
            else {
                return 0;
            }
        }
    }

    private static class ComparatorByTranscriptName implements Comparator<GencodeFuncotation> {

        public ComparatorByTranscriptName(){}

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // Compare ABC order by transcript name:
            if ( (a.getAnnotationTranscript() != null) && (b.getAnnotationTranscript() == null) ) {
                return -1;
            }
            else if ( (a.getAnnotationTranscript() == null ) && (b.getAnnotationTranscript() != null) ) {
                return 1;
            }
            // Need a default case in case all the comparison criteria are the same:
            else if ( (a.getAnnotationTranscript() == null ) && (b.getAnnotationTranscript() == null) ) {
                return -1;
            }
            else {
                return a.getAnnotationTranscript().compareTo(b.getAnnotationTranscript());
            }
        }
    }

    /**
     * Comparator class for Best Effect order.
     * Complex enough that a Lambda would be utter madness.
     */
    static class BestEffectGencodeFuncotationComparator implements Comparator<GencodeFuncotation> {

        private final Comparator<GencodeFuncotation> byUserTranscript;
        private final Comparator<GencodeFuncotation> byIgrStatus;
        private final Comparator<GencodeFuncotation> byVariantClassification;
        private final Comparator<GencodeFuncotation> byProteinCodingStatus;
        private final Comparator<GencodeFuncotation> byLocusLevel;
        private final Comparator<GencodeFuncotation> byApprisRank;
        private final Comparator<GencodeFuncotation> byTranscriptLength;
        private final Comparator<GencodeFuncotation> byTranscriptName;

        private final Comparator<GencodeFuncotation> chainedComparator;

        public BestEffectGencodeFuncotationComparator( final Set<String> userRequestedTranscripts ) {
            byUserTranscript = new ComparatorByUserTranscript(userRequestedTranscripts);
            byIgrStatus = new ComparatorByIgrStatus();
            byVariantClassification = new ComparatorByVariantClassification();
            byProteinCodingStatus = new ComparatorByProteinCodingStatus();
            byLocusLevel = new ComparatorByLocusLevel();
            byApprisRank = new ComparatorByApprisRank();
            byTranscriptLength = new ComparatorByTranscriptSequenceLength();
            byTranscriptName = new ComparatorByTranscriptName();

            chainedComparator = byUserTranscript
                    .thenComparing(byIgrStatus)
                    .thenComparing(byVariantClassification)
                    .thenComparing(byProteinCodingStatus)
                    .thenComparing(byLocusLevel)
                    .thenComparing(byApprisRank)
                    .thenComparing(byTranscriptLength)
                    .thenComparing(byTranscriptName);
        }

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // 0   - Choose the transcript that is on the custom list specified by the user
            // 1a) - Check to see if one is an IGR.
            // 1b) - Check highest variant classification.
            // 2a) - Is it protein coding?
            // 2b) - Check locus/curation levels.
            // 3)  - Check the appris annotation.
            // 4)  - Check transcript sequence length.
            // 5)  - Default to ABC order by transcript name.
            return chainedComparator.compare(a, b);
        }
    }

    /**
     * Comparator class for Canonical order.
     * Complex enough that a Lambda would be utter madness.
     */
    static class CanonicalGencodeFuncotationComparator implements Comparator<GencodeFuncotation> {

        private final Comparator<GencodeFuncotation> byUserTranscript;
        private final Comparator<GencodeFuncotation> byIgrStatus;
        private final Comparator<GencodeFuncotation> byVariantClassification;
        private final Comparator<GencodeFuncotation> byProteinCodingStatus;
        private final Comparator<GencodeFuncotation> byLocusLevel;
        private final Comparator<GencodeFuncotation> byApprisRank;
        private final Comparator<GencodeFuncotation> byTranscriptLength;
        private final Comparator<GencodeFuncotation> byTranscriptName;

        private final Comparator<GencodeFuncotation> chainedComparator;

        public CanonicalGencodeFuncotationComparator(final Set<String> userRequestedTranscripts ) {
            byUserTranscript = new ComparatorByUserTranscript(userRequestedTranscripts);
            byIgrStatus = new ComparatorByIgrStatus();
            byVariantClassification = new ComparatorByVariantClassification();
            byProteinCodingStatus = new ComparatorByProteinCodingStatus();
            byLocusLevel = new ComparatorByLocusLevel();
            byApprisRank = new ComparatorByApprisRank();
            byTranscriptLength = new ComparatorByTranscriptSequenceLength();
            byTranscriptName = new ComparatorByTranscriptName();

            chainedComparator = byUserTranscript
                    .thenComparing(byProteinCodingStatus)
                    .thenComparing(byLocusLevel)
                    .thenComparing(byApprisRank)
                    .thenComparing(byIgrStatus)
                    .thenComparing(byVariantClassification)
                    .thenComparing(byTranscriptLength)
                    .thenComparing(byTranscriptName);
        }

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {

            // 0)  - Choose the transcript that is on the custom list specified by the user:
            // 1a) - Is it protein coding?
            // 1b) - Check locus/curation levels.
            // 2)  - Check the appris annotation:
            // 3a) - Check to see if one is an IGR.
            // 3b) - Check highest variant classification:
            // 4)  - Check transcript sequence length:
            // 5)  - Default to ABC order by transcript name:
            return chainedComparator.compare(a, b);
        }
    }
}
