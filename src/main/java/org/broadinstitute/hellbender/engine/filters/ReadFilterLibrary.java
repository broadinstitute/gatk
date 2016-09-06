package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.Cigar;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Standard ReadFilters
 */
public final class ReadFilterLibrary {

    private ReadFilterLibrary(){ /*no instance*/ }

    /**
     * local classes for static read filters
     */
    public static class AllowAllReadsReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){return true;}}

    //Note: do not call getCigar to avoid creation of new Cigar objects
    public static class CigarIsSupportedReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! CigarUtils.containsNOperator(read.getCigarElements());}}

    public static class GoodCigarReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test (final GATKRead read) {
            return CigarUtils.isGood(read.getCigar());}}

    public static class HasMatchingBasesAndQualsReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getLength() == read.getBaseQualityCount();}}

    public static class HasReadGroupReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getReadGroup() != null;}}

    public static class MappedReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return !read.isUnmapped();}}

    public static class MappingQualityAvailableReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;}}

    public static class MappingQualityNotZeroReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getMappingQuality() != 0;}}

    /**
     * Reads that either have a mate that maps to the same contig, or don't have a mapped mate.
     */
    public static class MateOnSameContigOrNoMappedMateReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! read.isPaired() ||
                    read.mateIsUnmapped() ||
                    (! read.isUnmapped() && read.getContig().equals(read.getMateContig()));}}

    /**
     * Reads that have a mapped mate and both mate and read are on different same strands (ie the usual situation).
     */
    public static class MateDifferentStrandReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.isPaired() &&
                    ! read.isUnmapped() &&
                    ! read.mateIsUnmapped() &&
                    read.mateIsReverseStrand() != read.isReverseStrand();}}

    public static class NonZeroReferenceLengthAlignmentReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test (final GATKRead read) {
            return read.getCigarElements()
                    .stream()
                    .anyMatch(c -> c.getOperator().consumesReferenceBases() && c.getLength() > 0);
        }}

    public static class NotDuplicateReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! read.isDuplicate();}}

    public static class PassesVendorQualityCheckReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! read.failsVendorQualityCheck();}}

    public static class PrimaryAlignmentReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return ! read.isSecondaryAlignment();}}

    //Note: do not call getCigar to avoid creation of new Cigar objects
    public static class ReadLengthEqualsCigarLengthReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test (final GATKRead read) {
            return read.isUnmapped() ||
                    read.getLength() == Cigar.getReadLength(read.getCigarElements());}}

    public static class SeqIsStoredReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){
            return read.getLength() > 0;}}

    public static class ValidAlignmentStartReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read){return read.isUnmapped() || read.getStart() > 0;}}

    // Alignment doesn't align to a negative number of bases in the reference.
    //Note: to match gatk3 we must keep reads that align to zero bases in the reference.
    public static class ValidAlignmentEndReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final GATKRead read) {
            return read.isUnmapped() || (read.getEnd() - read.getStart() + 1) >= 0;}}

    /**
     * Static, stateless read filter instances
     */
    public static final AllowAllReadsReadFilter ALLOW_ALL_READS = new AllowAllReadsReadFilter();
    public static final CigarIsSupportedReadFilter CIGAR_IS_SUPPORTED = new CigarIsSupportedReadFilter();
    public static final GoodCigarReadFilter GOOD_CIGAR = new GoodCigarReadFilter();
    public static final HasReadGroupReadFilter HAS_READ_GROUP = new HasReadGroupReadFilter();
    public static final HasMatchingBasesAndQualsReadFilter HAS_MATCHING_BASES_AND_QUALS = new HasMatchingBasesAndQualsReadFilter();
    public static final MappedReadFilter MAPPED = new MappedReadFilter();
    public static final MappingQualityAvailableReadFilter MAPPING_QUALITY_AVAILABLE = new MappingQualityAvailableReadFilter();
    public static final MappingQualityNotZeroReadFilter MAPPING_QUALITY_NOT_ZERO   = new MappingQualityNotZeroReadFilter();
    public static final MateOnSameContigOrNoMappedMateReadFilter MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE = new MateOnSameContigOrNoMappedMateReadFilter();
    public static final MateDifferentStrandReadFilter MATE_DIFFERENT_STRAND = new MateDifferentStrandReadFilter();
    public static final NotDuplicateReadFilter NOT_DUPLICATE = new NotDuplicateReadFilter();
    public static final NonZeroReferenceLengthAlignmentReadFilter NON_ZERO_REFERENCE_LENGTH_ALIGNMENT = new NonZeroReferenceLengthAlignmentReadFilter();
    public static final PassesVendorQualityCheckReadFilter PASSES_VENDOR_QUALITY_CHECK = new PassesVendorQualityCheckReadFilter();
    public static final PrimaryAlignmentReadFilter PRIMARY_ALIGNMENT = new PrimaryAlignmentReadFilter();
    public static final ReadLengthEqualsCigarLengthReadFilter READLENGTH_EQUALS_CIGARLENGTH = new ReadLengthEqualsCigarLengthReadFilter();
    public static final SeqIsStoredReadFilter SEQ_IS_STORED = new SeqIsStoredReadFilter();
    public static final ValidAlignmentStartReadFilter VALID_ALIGNMENT_START = new ValidAlignmentStartReadFilter();
    public static final ValidAlignmentEndReadFilter VALID_ALIGNMENT_END = new ValidAlignmentEndReadFilter();

}
