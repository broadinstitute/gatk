package org.broadinstitute.hellbender.utils.smithwaterman;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.List;

/**
 * Pairwise discrete smith-waterman alignment
 *
 * ************************************************************************
 * ****                    IMPORTANT NOTE:                             ****
 * ****  This class assumes that all bytes come from UPPERCASED chars! ****
 * ************************************************************************
 */
public final class SWPairwiseAlignment {

    public static final SmithWatermanParameters STANDARD_NGS = new SmithWatermanParameters(25, -50, -110, -6, OverhangStrategy.SOFTCLIP);

    private final CoreSmithWatermanAlgorithm kernel;
    private final OverhangStrategy overhangStrategy;

    /**
     * Create a new SW pairwise aligner
     *
     * After creating the object the two sequences are aligned with an internal call to align(seq1, seq2)
     *
     * @param seq1 the first sequence we want to align
     * @param seq2 the second sequence we want to align
     * @param parameters the SW parameters to use
     */
    public SWPairwiseAlignment(final byte[] seq1, final byte[] seq2, final SmithWatermanParameters parameters) {
        this.kernel = new NativeSSECoreSmithWatermanAlgorithm();//new JavaCoreSmithWatermanImplementation();
        final boolean ok = kernel.initialize(parameters);
        overhangStrategy = parameters.overhangStrategy;
        kernel.align(seq1, seq2);
    }

    public Cigar getCigar() { return kernel.getCigar() ; }

    public int getAlignmentStart2wrt1() { return kernel.getAlignmentStart2wrt1(); }

    @VisibleForTesting
    void printAlignment(final byte[] ref, final byte[] read) {
        printAlignment(ref,read,100);
    }

    @VisibleForTesting
    public void printAlignment(final byte[] ref, final byte[] read, final int width) {
        final StringBuilder bread = new StringBuilder();
        final StringBuilder bref = new StringBuilder();
        final StringBuilder match = new StringBuilder();

        int i = 0;
        int j = 0;

        final int offset = getAlignmentStart2wrt1();

        Cigar cigar = getCigar();

        if ( overhangStrategy != OverhangStrategy.SOFTCLIP ) {

            // we need to go through all the hassle below only if we do not do softclipping;
            // otherwise offset is never negative
            if ( offset < 0 ) {
                for (  ; j < (-offset) ; j++ ) {
                    bread.append((char)read[j]);
                    bref.append(' ');
                    match.append(' ');
                }
                // at negative offsets, our cigar's first element carries overhanging bases
                // that we have just printed above. Tweak the first element to
                // exclude those bases. Here we create a new list of cigar elements, so the original
                // list/original cigar are unchanged (they are unmodifiable anyway!)

                final List<CigarElement> tweaked = new ArrayList<>();
                tweaked.addAll(cigar.getCigarElements());
                tweaked.set(0,new CigarElement(cigar.getCigarElement(0).getLength()+offset,
                        cigar.getCigarElement(0).getOperator()));
                cigar = new Cigar(tweaked);
            }
        }

        if ( offset > 0 ) { // note: the way this implementation works, cigar will ever start from S *only* if read starts before the ref, i.e. offset = 0
            for (  ; i < getAlignmentStart2wrt1() ; i++ ) {
                bref.append((char)ref[i]);
                bread.append(' ');
                match.append(' ');
            }
        }
        
        for ( final CigarElement e : cigar.getCigarElements() ) {
            switch (e.getOperator()) {
                case M :
                    for ( int z = 0 ; z < e.getLength() ; z++, i++, j++  ) {
                        bref.append((i<ref.length)?(char)ref[i]:' ');
                        bread.append((j < read.length)?(char)read[j]:' ');
                        match.append( ( i<ref.length && j < read.length ) ? (ref[i] == read[j] ? '.':'*' ) : ' ' );
                    }
                    break;
                case I :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append('-');
                        bread.append((char)read[j]);
                        match.append('I');
                    }
                    break;
                case S :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append(' ');
                        bread.append((char)read[j]);
                        match.append('S');
                    }
                    break;
                case D:
                    for ( int z = 0 ; z < e.getLength(); z++ , i++ ) {
                        bref.append((char)ref[i]);
                        bread.append('-');
                        match.append('D');
                    }
                    break;
                default:
                    throw new GATKException("Unexpected Cigar element:" + e.getOperator());
            }
        }
        for ( ; i < ref.length; i++ ) bref.append((char)ref[i]);
        for ( ; j < read.length; j++ ) bread.append((char)read[j]);

        int pos = 0 ;
        final int maxlength = Math.max(match.length(), Math.max(bread.length(), bref.length()));
        while ( pos < maxlength ) {
            print_cautiously(match,pos,width);
            print_cautiously(bread,pos,width);
            print_cautiously(bref,pos,width);
            System.out.println();
            pos += width;
        }
    }

    /** String builder's substring is extremely stupid: instead of trimming and/or returning an empty
     * string when one end/both ends of the interval are out of range, it crashes with an
     * exception. This utility function simply prints the substring if the interval is within the index range
     * or trims accordingly if it is not.
     * @param s
     * @param start
     * @param width
     */
    private static void print_cautiously(final StringBuilder s, final int start, final int width) {
        if ( start >= s.length() ) {
            System.out.println();
            return;
        }
        final int end = Math.min(start + width, s.length());
        System.out.println(s.substring(start,end));
    }
}
