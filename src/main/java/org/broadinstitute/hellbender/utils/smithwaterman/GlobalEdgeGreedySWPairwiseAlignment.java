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

package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.sam.AlignmentUtils;

import java.util.*;

/**
 * Pairwise discrete Smith-Waterman alignment with an edge greedy implementation
 *
 * ************************************************************************
 * ****                    IMPORTANT NOTE:                             ****
 * ****  This class assumes that all bytes come from UPPERCASED chars! ****
 * ************************************************************************
 *
 * User: ebanks
 */
public final class GlobalEdgeGreedySWPairwiseAlignment extends SWPairwiseAlignment {

    private final static boolean DEBUG_MODE = false;

    /**
     * Create a new greedy SW pairwise aligner
     *
     * @param reference the reference sequence we want to align
     * @param alternate the alternate sequence we want to align
     * @param parameters the SW parameters to use
     */
    public GlobalEdgeGreedySWPairwiseAlignment(final byte[] reference, final byte[] alternate, final Parameters parameters) {
        super(reference, alternate, parameters);
    }

    /**
     * Create a new SW pairwise aligner
     *
     * After creating the object the two sequences are aligned with an internal call to align(seq1, seq2)
     *
     * @param reference the reference sequence we want to align
     * @param alternate the alternate sequence we want to align
     * @param namedParameters the named parameter set to get our parameters from
     */
    public GlobalEdgeGreedySWPairwiseAlignment(final byte[] reference, final byte[] alternate, final SWParameterSet namedParameters) {
        this(reference, alternate, namedParameters.parameters);
    }

    /**
     * @see #GlobalEdgeGreedySWPairwiseAlignment(byte[], byte[], SWParameterSet) with original default parameters
     */
    public GlobalEdgeGreedySWPairwiseAlignment(byte[] reference, byte[] alternate) {
        this(reference, alternate, SWParameterSet.ORIGINAL_DEFAULT);
    }

    /**
     * Aligns the alternate sequence to the reference sequence
     *
     * @param reference  ref sequence
     * @param alternate  alt sequence
     */
    @Override
    protected void align(final byte[] reference, final byte[] alternate) {
        if ( reference == null || reference.length == 0 )
            throw new IllegalArgumentException("Non-null, non-empty reference sequences are required for the Smith-Waterman calculation");
        if ( alternate == null || alternate.length == 0 )
            throw new IllegalArgumentException("Non-null, non-empty alternate sequences are required for the Smith-Waterman calculation");

        final int forwardEdgeMatch = Utils.longestCommonPrefix(reference, alternate, Integer.MAX_VALUE);

        // edge case: one sequence is a strict prefix of the other
        if ( forwardEdgeMatch == reference.length || forwardEdgeMatch == alternate.length ) {
            alignmentResult = new SWPairwiseAlignmentResult(makeCigarForStrictPrefixAndSuffix(reference, alternate, forwardEdgeMatch, 0), 0);
            return;
        }

        int reverseEdgeMatch = Utils.longestCommonSuffix(reference, alternate, Integer.MAX_VALUE);

        // edge case: one sequence is a strict suffix of the other
        if ( reverseEdgeMatch == reference.length || reverseEdgeMatch == alternate.length ) {
            alignmentResult = new SWPairwiseAlignmentResult(makeCigarForStrictPrefixAndSuffix(reference, alternate, 0, reverseEdgeMatch), 0);
            return;
        }

        final int sizeOfRefToAlign = reference.length - forwardEdgeMatch - reverseEdgeMatch;
        final int sizeOfAltToAlign = alternate.length - forwardEdgeMatch - reverseEdgeMatch;

        // edge case: one sequence is a strict subset of the other accounting for both prefix and suffix
        final int minSizeToAlign = Math.min(sizeOfRefToAlign, sizeOfAltToAlign);
        if ( minSizeToAlign < 0 )
            reverseEdgeMatch += minSizeToAlign;
        if ( sizeOfRefToAlign <= 0 || sizeOfAltToAlign <= 0 ) {
            alignmentResult = new SWPairwiseAlignmentResult(makeCigarForStrictPrefixAndSuffix(reference, alternate, forwardEdgeMatch, reverseEdgeMatch), 0);
            return;
        }

        final byte[] refToAlign = Utils.trimArray(reference, forwardEdgeMatch, reverseEdgeMatch);
        final byte[] altToAlign = Utils.trimArray(alternate, forwardEdgeMatch, reverseEdgeMatch);

        final int[][] sw = new int[(sizeOfRefToAlign+1)][(sizeOfAltToAlign+1)];
        if ( keepScoringMatrix ) SW = sw;
        final int[][] btrack = new int[(sizeOfRefToAlign+1)][(sizeOfAltToAlign+1)];

        calculateMatrix(refToAlign, altToAlign, sw, btrack, OVERHANG_STRATEGY.INDEL);

        if ( DEBUG_MODE ) {
            System.out.println(new String(refToAlign) + " vs. " + new String(altToAlign));
            debugMatrix(sw);
            System.out.println("----");
            debugMatrix(btrack);
            System.out.println();
        }

        alignmentResult = calculateCigar(forwardEdgeMatch, reverseEdgeMatch, sw, btrack);
    }


    private void debugMatrix(final int[][] matrix) {
        for ( int i = 0; i < matrix.length; i++ ) {
            int [] cur = matrix[i];
            for ( int j = 0; j < cur.length; j++ )
                System.out.print(cur[j] + " ");
            System.out.println();
        }
    }

    /**
     * Creates a CIGAR for the case where the prefix/suffix match combination encompasses an entire sequence
     *
     * @param reference            the reference sequence
     * @param alternate            the alternate sequence
     * @param matchingPrefix       the prefix match size
     * @param matchingSuffix       the suffix match size
     * @return non-null CIGAR
     */
    private Cigar makeCigarForStrictPrefixAndSuffix(final byte[] reference, final byte[] alternate, final int matchingPrefix, final int matchingSuffix) {

        final List<CigarElement> result = new ArrayList<CigarElement>();

        // edge case: no D or I element
        if ( reference.length == alternate.length ) {
            result.add(makeElement(State.MATCH, matchingPrefix + matchingSuffix));
        } else {
            // add the first M element
            if ( matchingPrefix > 0 )
                result.add(makeElement(State.MATCH, matchingPrefix));

            // add the D or I element
            if ( alternate.length > reference.length )
                result.add(makeElement(State.INSERTION, alternate.length - reference.length));
            else // if ( reference.length > alternate.length )
                result.add(makeElement(State.DELETION, reference.length - alternate.length));

            // add the last M element
            if ( matchingSuffix > 0 )
                result.add(makeElement(State.MATCH, matchingSuffix));
        }

        return new Cigar(result);
    }

    /**
     * Calculates the CIGAR for the alignment from the back track matrix
     *
     * @param matchingPrefix       the prefix match size
     * @param matchingSuffix       the suffix match size
     * @param sw                   the Smith-Waterman matrix to use
     * @param btrack               the back track matrix to use
     * @return non-null SWPairwiseAlignmentResult object
     */
    protected SWPairwiseAlignmentResult calculateCigar(final int matchingPrefix, final int matchingSuffix,
                                                       final int[][] sw, final int[][] btrack) {

        final SWPairwiseAlignmentResult SW_result = calculateCigar(sw, btrack, OVERHANG_STRATEGY.INDEL);

        final LinkedList<CigarElement> lce = new LinkedList<CigarElement>(SW_result.cigar.getCigarElements());
        if ( matchingPrefix > 0 )
            lce.addFirst(makeElement(State.MATCH, matchingPrefix));
        if ( matchingSuffix > 0 )
            lce.addLast(makeElement(State.MATCH, matchingSuffix));

        return new SWPairwiseAlignmentResult(AlignmentUtils.consolidateCigar(new Cigar(lce)), 0);
    }
}