package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Tests whether a read is &quot;well-formed&quot; -- that is, is free of major internal inconsistencies and issues that could lead
 * to errors downstream. If a read passes this filter, the rest of the engine should be able to process it without
 * blowing up.
 *
 * <p><b>Well-formed reads definition</b></p>
 * <ul>
 *     <li><b>Alignment coordinates:</b> start larger than 0 and end after the start position.</li>
 *     <li><b>Alignment agrees with header:</B> contig exists and start is within its range.</li>
 *     <li><b>Read Group and Sequence are present</b></li>
 *     <li><b>Consistent read length:</b> bases match in length with the qualities and the CIGAR string.</b></li>
 *     <li><b>Do not contain skipped regions:</b> represented by the 'N' operator in the CIGAR string.</li>
 * </ul>
 *
 * @see ReadFilterLibrary.ValidAlignmentStartReadFilter
 * @see ReadFilterLibrary.ValidAlignmentEndReadFilter
 * @see AlignmentAgreesWithHeaderReadFilter
 * @see ReadFilterLibrary.HasReadGroupReadFilter
 * @see ReadFilterLibrary.MatchingBasesAndQualsReadFilter
 * @see ReadFilterLibrary.ReadLengthEqualsCigarLengthReadFilter
 * @see ReadFilterLibrary.SeqIsStoredReadFilter
 * @see ReadFilterLibrary.CigarContainsNoNOperator
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads that are well-formed",
        extraDocs = {
                ReadFilterLibrary.ValidAlignmentStartReadFilter.class,
                ReadFilterLibrary.ValidAlignmentEndReadFilter.class,
                AlignmentAgreesWithHeaderReadFilter.class,
                ReadFilterLibrary.HasReadGroupReadFilter.class,
                ReadFilterLibrary.MatchingBasesAndQualsReadFilter.class,
                ReadFilterLibrary.ReadLengthEqualsCigarLengthReadFilter.class,
                ReadFilterLibrary.SeqIsStoredReadFilter.class,
                ReadFilterLibrary.CigarContainsNoNOperator.class
        }
)
public final class WellformedReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1l;

    private ReadFilter wellFormedFilter = null;

    // Command line parser requires a no-arg constructor
    public  WellformedReadFilter() {
    }

    @Override
    public void setHeader(SAMFileHeader header) {
        super.setHeader(header);
        createFilter();
    }

    public WellformedReadFilter(final SAMFileHeader header) {
        setHeader(header);
    }

    private void createFilter() {
        final AlignmentAgreesWithHeaderReadFilter alignmentAgreesWithHeader = new AlignmentAgreesWithHeaderReadFilter(samHeader);

        wellFormedFilter = ReadFilterLibrary.VALID_ALIGNMENT_START
                .and(ReadFilterLibrary.VALID_ALIGNMENT_END)
                .and(alignmentAgreesWithHeader)
                .and(ReadFilterLibrary.HAS_READ_GROUP)
                .and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                .and(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH)
                .and(ReadFilterLibrary.SEQ_IS_STORED)
                .and(ReadFilterLibrary.CIGAR_CONTAINS_NO_N_OPERATOR);
    }

    @Override
    public boolean test(final GATKRead read ) {
        return wellFormedFilter.test(read);
    }
}
