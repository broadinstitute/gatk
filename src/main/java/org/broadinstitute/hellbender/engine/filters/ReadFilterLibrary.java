package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.Cigar;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Standard ReadFilters
 */
public final class ReadFilterLibrary {

    private ReadFilterLibrary(){ /*no instance*/ }

    /** Do not filter out any read. */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Do not filter out any read")
    public static class AllowAllReadsReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){return true;}}

    /** Filter out reads containing skipped region from the reference (CIGAR strings with 'N' operator). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads with CIGAR containing N operator")
    //Note: do not call getCigar to avoid creation of new Cigar objects
    public static class CigarContainsNoNOperator extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! CigarUtils.containsNOperator(read.getCigarElements());}}

    /** Keep only reads that are first of pair (0x1 and 0x40). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads that are first of pair")
    public static class FirstOfPairReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test (final GATKRead read) {
            return read.isFirstOfPair();}}

    /**
     * Keep only reads containing good CIGAR strings.
     *
     * <p>Good CIGAR strings have the following properties:</p>
     *
     * <ul>
     *     <li>Valid according to the <a href="http://samtools.github.io/hts-specs/SAMv1.pdf">SAM specifications.</a></li>
     *     <li>Does not start or end with deletions (with or without preceding clips).</li>
     *     <li>Does not have consecutive deletion/insertion operators.</li>
     * </ul>
     */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads containing good CIGAR string")
    public static class GoodCigarReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test (final GATKRead read) {
            return CigarUtils.isGood(read.getCigar());}}

    /** Filter out reads with fragment length (insert size) different from zero. */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads with fragment length different from zero")
    public static class NonZeroFragmentLengthReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getFragmentLength() != 0;}}

    /** Filter out reads where the bases and qualities do not match in length. */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads where the bases and qualities do not match")
    public static class MatchingBasesAndQualsReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getLength() == read.getBaseQualityCount();}}

    /** Filter out reads without the SAM record RG (Read Group) tag. */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads without Read Group")
    public static class HasReadGroupReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getReadGroup() != null;}}

    /**
     * Filter out unmapped reads.
     *
     * <p>Unmapped reads are defined by three criteria:</p>
     *
     * <ul>
     *     <li>SAM flag value 0x4</li>
     *     <li>An asterisk for reference name or RNAME (column 3 of a SAM record)</li>
     *     <li>A zero value for leftmost mapping position for POS (column 4 of SAM record)</li>
     * </ul>
     */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out unmapped reads")
    public static class MappedReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return !read.isUnmapped();}}

    /** Filter out reads without available mapping quality (MAPQ=255). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads without available mapping quality")
    public static class MappingQualityAvailableReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;}}

    /** Filter out reads with mapping quality equal to zero. */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads with mapping quality equal to zero")
    public static class MappingQualityNotZeroReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return read.getMappingQuality() != 0;}}

    /**
     * Keep only reads that have a mate that maps to the same contig (RNEXT is "="), is single ended (not 0x1) or has an unmapped mate (0x8).
     *
     * <p>See MappedReadFilter for criteria defining an unmapped read.</p>
     */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads whose mate maps to the same contig or is unmapped", extraDocs = MappedReadFilter.class)
    public static class MateOnSameContigOrNoMappedMateReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! read.isPaired() ||
                    read.mateIsUnmapped() ||
                    (! read.isUnmapped() && read.getContig().equals(read.getMateContig()));}}

    /**
     * For paired reads (0x1), keep only reads that are mapped, have a mate that is mapped (read is not 0x8), and both
     * the read and its mate are on different strands (when read is 0x20, it is not 0x10), as is the typical case.
     *
     * <p>See MappedReadFilter for criteria defining an mapped read.</p>
     */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads with mates mapped on the different strand", extraDocs = MappedReadFilter.class)
    public static class MateDifferentStrandReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.isPaired() &&
                    ! read.isUnmapped() &&
                    ! read.mateIsUnmapped() &&
                    read.mateIsReverseStrand() != read.isReverseStrand();}}

    /** Filter out reads that do not align to the reference. Filter interprets each of the CIGAR operators M, D, N, = and X as alignment. */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads that do not align to the reference")
    public static class NonZeroReferenceLengthAlignmentReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test (final GATKRead read) {
            return read.getCigarElements()
                    .stream()
                    .anyMatch(c -> c.getOperator().consumesReferenceBases() && c.getLength() > 0);
        }}

    /** Filter out reads marked as duplicate (0x400). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads marked as duplicate")
    public static class NotDuplicateReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! read.isDuplicate();}}

    /** Filter out reads representing secondary alignments (0x100). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads representing secondary alignments")
    public static class NotSecondaryAlignmentReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return !read.isSecondaryAlignment();}}

    /** Filter out reads representing supplementary alignments (0x800). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads representing supplementary alignments")
    public static class NotSupplementaryAlignmentReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return !read.isSupplementaryAlignment();}}

    /** Filter out unpaired reads (not 0x1). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out unpaired reads")
    public static class PairedReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return read.isPaired();}}

    /** Filter out reads failing platform/vendor quality checks (0x200). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads failing platfor/vendor quality checks")
    public static class PassesVendorQualityCheckReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! read.failsVendorQualityCheck();}}

    /** Keep only paired reads that are properly paired (0x1 and 0x2). Removes single ended reads. */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads that are properly paired")
    public static class ProperlyPairedReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return read.isProperlyPaired();}}

    /**
     * Keep only reads representing primary alignments (those that satisfy both the NotSecondaryAlignment and
     * NotSupplementaryAlignment filters, or in terms of SAM flag values, must have neither of the 0x100 or
     * 0x800 flags set).
     *
     * <p>Note that this filter represents a stronger criteria for "primary alignment" than the
     * SAM flag 0x100 (representing ""not primary alignment" in some contexts).</p>
     *
     * <p>For example, a read that has only the supplementary flag (0x800) set, but not the secondary (0x100)
     * flag will be filtered out from processing by the PrimaryLineReadFilter, but would NOT be filtered out by
     * other software that uses the looser notion of "not primary" that only depends on the "secondary" flag being set.</p>
     */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY)
    public static class PrimaryLineReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return ! read.isSecondaryAlignment() && ! read.isSupplementaryAlignment();}}

    /**
     * Filter out reads where the read and CIGAR do not match in length.
     *
     * <p>Note: unmapped reads pass this filter. See MappedReadFilter for criteria defining an unmapped read.</p>
     */
    //Note: do not call getCigar to avoid creation of new Cigar objects
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads  where the read and CIGAR do not match in length", extraDocs = MappedReadFilter.class)
    public static class ReadLengthEqualsCigarLengthReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test (final GATKRead read) {
            return read.isUnmapped() ||
                    read.getLength() == Cigar.getReadLength(read.getCigarElements());}}

    /** Keep only paired reads (0x1) that are second of pair (0x80). */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only paired reads that are second of pair")
    public static class SecondOfPairReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test (final GATKRead read) {
            return read.isSecondOfPair();}}

    /** Keep only reads with sequenced bases. */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads with sequenced bases")
    public static class SeqIsStoredReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getLength() > 0;}}

    /**
     * Keep only reads with a valid alignment start (POS larger than 0) or is unmapped.
     *
     * <p>See MappedReadFilter for criteria defining an unmapped read.</p>
     */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads with a valid alignment start", extraDocs = MappedReadFilter.class)
    public static class ValidAlignmentStartReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.isUnmapped() || read.getStart() > 0;}}

    /**
     * Keep only reads where the read end corresponds to a proper alignment -- that is, the read ends after the start
     * (non-negative number of bases in the reference).
     *
     * <p>This is calculated as:</p>
     *
     * <p>
     * <code>
     * end - start + 1 &ge; 0<br>
     * where<br>
     *  start = 1-based inclusive leftmost position of the clipped sequence (0 if no position)<br>
     *  end = 1-based inclusive rightmost position of the clipped sequence (0 if unmapped)<br>
     * </code>
     * </p>
     *
     * <p>Note: keep also unmapped reads (align to zero bases in the reference). See MappedReadFilter for criteria defining an unmapped read.</p>
     */
    @DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads where the read end is properly aligned", extraDocs = MappedReadFilter.class)
    public static class ValidAlignmentEndReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return read.isUnmapped() || (read.getEnd() - read.getStart() + 1) >= 0;}}

    /**
     * Static, stateless read filter instances
     */
    public static final AllowAllReadsReadFilter ALLOW_ALL_READS = new AllowAllReadsReadFilter();
    public static final CigarContainsNoNOperator CIGAR_CONTAINS_NO_N_OPERATOR = new CigarContainsNoNOperator();
    public static final FirstOfPairReadFilter FIRST_OF_PAIR = new FirstOfPairReadFilter();
    public static final GoodCigarReadFilter GOOD_CIGAR = new GoodCigarReadFilter();
    public static final HasReadGroupReadFilter HAS_READ_GROUP = new HasReadGroupReadFilter();
    public static final MappedReadFilter MAPPED = new MappedReadFilter();
    public static final MappingQualityAvailableReadFilter MAPPING_QUALITY_AVAILABLE = new MappingQualityAvailableReadFilter();
    public static final MappingQualityNotZeroReadFilter MAPPING_QUALITY_NOT_ZERO   = new MappingQualityNotZeroReadFilter();
    public static final MatchingBasesAndQualsReadFilter HAS_MATCHING_BASES_AND_QUALS = new MatchingBasesAndQualsReadFilter();
    public static final MateOnSameContigOrNoMappedMateReadFilter MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE = new MateOnSameContigOrNoMappedMateReadFilter();
    public static final MateDifferentStrandReadFilter MATE_DIFFERENT_STRAND = new MateDifferentStrandReadFilter();
    public static final NonZeroReferenceLengthAlignmentReadFilter NON_ZERO_REFERENCE_LENGTH_ALIGNMENT = new NonZeroReferenceLengthAlignmentReadFilter();
    public static final NonZeroFragmentLengthReadFilter NONZERO_FRAGMENT_LENGTH_READ_FILTER = new NonZeroFragmentLengthReadFilter();
    public static final NotDuplicateReadFilter NOT_DUPLICATE = new NotDuplicateReadFilter();
    public static final NotSecondaryAlignmentReadFilter NOT_SECONDARY_ALIGNMENT = new NotSecondaryAlignmentReadFilter();
    public static final NotSupplementaryAlignmentReadFilter NOT_SUPPLEMENTARY_ALIGNMENT = new NotSupplementaryAlignmentReadFilter();
    public static final PairedReadFilter PAIRED = new PairedReadFilter();
    public static final ProperlyPairedReadFilter PROPERLY_PAIRED = new ProperlyPairedReadFilter();
    public static final PassesVendorQualityCheckReadFilter PASSES_VENDOR_QUALITY_CHECK = new PassesVendorQualityCheckReadFilter();
    public static final PrimaryLineReadFilter PRIMARY_LINE = new PrimaryLineReadFilter();
    public static final ReadLengthEqualsCigarLengthReadFilter READLENGTH_EQUALS_CIGARLENGTH = new ReadLengthEqualsCigarLengthReadFilter();
    public static final SecondOfPairReadFilter SECOND_OF_PAIR = new SecondOfPairReadFilter();
    public static final SeqIsStoredReadFilter SEQ_IS_STORED = new SeqIsStoredReadFilter();
    public static final ValidAlignmentStartReadFilter VALID_ALIGNMENT_START = new ValidAlignmentStartReadFilter();
    public static final ValidAlignmentEndReadFilter VALID_ALIGNMENT_END = new ValidAlignmentEndReadFilter();

}
