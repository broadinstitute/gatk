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

/**
 * Defines what it means to be a wellformed read..
 */
public final class WellformedReads  {

    public static boolean isWellformed(final SAMRecord read) {
        return hasValidAlignmentStart(read)
                && hasValidAlignmentEnd(read)
                && alignmentAgreesWithHeader(read.getHeader(), read)
                && hasReadGroup(read)
                && hasMatchingBasesAndQuals(read)
                && cigarAgreesWithAlignment(read)
                && seqIsStored(read)
                && cigarIsSupported(read);
    }

    private static boolean hasReadGroup(final SAMRecord read) {
        return read.getReadGroup() != null;
    }

    /**
     * Check for the case in which the alignment start is inconsistent with the read unmapped flag.
     * @param read The read to validate.
     * @return true if read start is valid, false otherwise.
     */
    private static boolean hasValidAlignmentStart(final SAMRecord read) {
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
    private static boolean hasValidAlignmentEnd(final SAMRecord read) {
        // Alignment aligns to negative number of bases in the reference.
        return !(!read.getReadUnmappedFlag() && read.getAlignmentEnd() != -1 && (read.getAlignmentEnd() - read.getAlignmentStart() + 1) < 0);
    }

    /**
     * Check to ensure that the alignment makes sense based on the contents of the header.
     * @param header The SAM file header.
     * @param read The read to verify.
     * @return true if alignment agrees with header, false otherwise.
     */
    private static boolean alignmentAgreesWithHeader(final SAMFileHeader header, final SAMRecord read) {
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
    private static boolean cigarAgreesWithAlignment(final SAMRecord read) {
        // Read has a valid alignment start, but the CIGAR string is empty
        return read.getReadUnmappedFlag() || read.getAlignmentStart() == -1 || read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START || read.getAlignmentBlocks().size() >= 0;
    }

    /**
     * Check for unsupported CIGAR operators.
     * Currently the N operator is not supported.
     * @param read The read to validate.
     * @return <code>true</code> if the read CIGAR operations are
     * fully supported, otherwise <code>false</code>.
     */
    private static boolean cigarIsSupported(final SAMRecord read) {
        return !containsNOperator(read);
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
    private static boolean hasMatchingBasesAndQuals(final SAMRecord read) {
        return read.getReadLength() == read.getBaseQualities().length;
    }

    /**
     * Check if the read has its base sequence stored
     * @param read the read to validate
     * @return true if the sequence is stored and false otherwise ("*" in the SEQ field).
     */
    protected static boolean seqIsStored(final SAMRecord read) {
        return read.getReadBases() != SAMRecord.NULL_SEQUENCE ;
    }
}
