/*
* Copyright (c) 2012 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Filter out malformed reads.
 */
public class MalformedReadFilter implements ReadFilter {

    /**
     * Filter out reads with CIGAR containing the N operator, instead of stop processing and report an error.
     */
    private final boolean filterReadsWithNCigar;

    /**
     * If a read has mismatching number of bases and base qualities, filter out the read instead of blowing up.
     */
    private final boolean filterMismatchingBaseAndQuals;

    /**
     * If a read has no stored bases (i.e. a '*'), filter out the read instead of blowing up.
     */
    private final boolean filterBasesNotStored;

    /**
     * Indicates whether to blow up on reads with a N in them.
     */
    private final boolean allowNCigars;

    public MalformedReadFilter(boolean filterReadsWithNCigar, boolean filterMismatchingBaseAndQuals, boolean filterBasesNotStored, boolean allowNCigars){
        this.filterReadsWithNCigar = filterReadsWithNCigar;
        this.filterMismatchingBaseAndQuals = filterMismatchingBaseAndQuals;
        this.filterBasesNotStored=filterBasesNotStored;
        this.allowNCigars=allowNCigars;
    }

    public MalformedReadFilter() {
       this(false, false, false, false);
    }

    @Override
    public boolean test(final SAMRecord read) {
        // slowly changing the behavior to blow up first and filtering out if a parameter is explicitly provided
        return  !checkInvalidAlignmentStart(read) ||
                !checkInvalidAlignmentEnd(read) ||
                !checkAlignmentDisagreesWithHeader(read.getHeader(),read) ||
                !checkHasReadGroup(read) ||
                !checkMismatchingBasesAndQuals(read, filterMismatchingBaseAndQuals) ||
                !checkCigarDisagreesWithAlignment(read) ||
                !checkSeqStored(read, filterBasesNotStored) ||
                !checkCigarIsSupported(read,filterReadsWithNCigar,allowNCigars);
    }

    private static boolean checkHasReadGroup(final SAMRecord read) {
        if ( read.getReadGroup() == null ) {
            // there are 2 possibilities: either the RG tag is missing or it is not defined in the header
            final String rgID = (String)read.getAttribute(SAMTagUtil.getSingleton().RG);
            if ( rgID == null )
                throw new UserException.ReadMissingReadGroup(read);
            throw new UserException.ReadHasUndefinedReadGroup(read, rgID);
        }
        return true;
    }

    /**
     * Check for the case in which the alignment start is inconsistent with the read unmapped flag.
     * @param read The read to validate.
     * @return true if read start is valid, false otherwise.
     */
    private static boolean checkInvalidAlignmentStart(final SAMRecord read ) {
        // read is not flagged as 'unmapped', but alignment start is NO_ALIGNMENT_START
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START )
            return false;
        // Read is not flagged as 'unmapped', but alignment start is -1
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() == -1 )
            return false;
        return true;
    }

    /**
     * Check for invalid end of alignments.
     * @param read The read to validate.
     * @return true if read end is valid, false otherwise.
     */
    private static boolean checkInvalidAlignmentEnd(final SAMRecord read ) {
        // Alignment aligns to negative number of bases in the reference.
        if( !read.getReadUnmappedFlag() && read.getAlignmentEnd() != -1 && (read.getAlignmentEnd()-read.getAlignmentStart()+1)<0 )
            return false;
        return true;
    }

    /**
     * Check to ensure that the alignment makes sense based on the contents of the header.
     * @param header The SAM file header.
     * @param read The read to verify.
     * @return true if alignment agrees with header, false otherwise.
     */
    private static boolean checkAlignmentDisagreesWithHeader(final SAMFileHeader header, final SAMRecord read ) {
        // Read is aligned to nonexistent contig
        if( read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START )
            return false;
        final SAMSequenceRecord contigHeader = header.getSequence( read.getReferenceIndex() );
        // Read is aligned to a point after the end of the contig
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() > contigHeader.getSequenceLength() )
            return false;
        return true;
    }

    /**
     * Check for inconsistencies between the cigar string and the
     * @param read The read to validate.
     * @return true if cigar agrees with alignment, false otherwise.
     */
    private static boolean checkCigarDisagreesWithAlignment(final SAMRecord read) {
        // Read has a valid alignment start, but the CIGAR string is empty
        if( !read.getReadUnmappedFlag() &&
                read.getAlignmentStart() != -1 &&
                read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START &&
                read.getAlignmentBlocks().size() < 0 )
            return false;
        return true;
    }

    /**
     * Check for unsupported CIGAR operators.
     * Currently the N operator is not supported.
     * @param read The read to validate.
     * @param filterReadsWithNCigar whether the offending read should just
     *                              be silently filtered or not.
     * @param allowNCigars whether reads that contain N operators in their CIGARs
     *                     can be processed or an exception should be thrown instead.
     * @throws UserException.UnsupportedCigarOperatorException
     *   if {@link #filterReadsWithNCigar} is <code>false</code> and
     *   the input read has some unsupported operation.
     * @return <code>true</code> if the read CIGAR operations are
     * fully supported, otherwise <code>false</code>, as long as
     * no exception has been thrown.
     */
    private static boolean checkCigarIsSupported(final SAMRecord read, final boolean filterReadsWithNCigar, final boolean allowNCigars) {
        if( containsNOperator(read)) {
            if (! filterReadsWithNCigar && !allowNCigars) {
                throw new UserException.UnsupportedCigarOperatorException(
                        CigarOperator.N,read,
                        "Perhaps you are"
                                + " trying to use RNA-Seq data?"
                                + " While we are currently actively working to"
                                + " support this data type unfortunately the"
                                + " GATK cannot be used with this data in its"
                                + " current form."); //TODO re-enable this in hellbender: https://github.com/broadinstitute/hellbender/issues/186
            }
            return ! filterReadsWithNCigar;
        }
        return true;
    }

    private static boolean containsNOperator(final SAMRecord read) {
        final Cigar cigar = read.getCigar();
        if (cigar == null)   {
            return false;
        }
        for (final CigarElement ce : cigar.getCigarElements()) {
            if (ce.getOperator() == CigarOperator.N) {
                return true;
            }
        }
        return false;
    }

    /**
     * Check if the read has the same number of bases and base qualities
     * @param read the read to validate
     * @return true if they have the same number. False otherwise.
     */
    private static boolean checkMismatchingBasesAndQuals(final SAMRecord read, final boolean filterMismatchingBaseAndQuals) {
        final boolean result;
        if (read.getReadLength() == read.getBaseQualities().length)
            result = true;
        else if (filterMismatchingBaseAndQuals)
            result = false;
        else
            throw new UserException.MalformedBAM(read,
                    String.format("BAM file has a read with mismatching number of bases and base qualities. Offender: %s [%d bases] [%d quals].",
                            read.getReadName(), read.getReadLength(), read.getBaseQualities().length));

        return result;
    }

    /**
     * Check if the read has its base sequence stored
     * @param read the read to validate
     * @return true if the sequence is stored and false otherwise ("*" in the SEQ field).
     */
    protected static boolean checkSeqStored(final SAMRecord read, final boolean filterBasesNotStored) {

        if ( read.getReadBases() != SAMRecord.NULL_SEQUENCE )
            return true;

        if ( filterBasesNotStored )
            return false;

        throw new UserException.MalformedBAM(read, String.format("the BAM file has a read with no stored bases (i.e. it uses '*') which is not supported in the GATK; see the --filter_bases_not_stored argument. Offender: %s", read.getReadName()));
    }
}
