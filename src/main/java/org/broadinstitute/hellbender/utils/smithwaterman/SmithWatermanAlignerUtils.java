package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.List;

/**
 * Utils for printing alignments
 */
public class SmithWatermanAlignerUtils {
    public static void printAlignment(final byte[] ref, final byte[] read, final SmithWatermanAlignment alignment, final SWOverhangStrategy overhangStrategy) {
        final StringBuilder bread = new StringBuilder();
        final StringBuilder bref = new StringBuilder();
        final StringBuilder match = new StringBuilder();

        int i = 0;
        int j = 0;

        final int offset = alignment.getAlignmentOffset();

        Cigar cigar = alignment.getCigar();

        if ( overhangStrategy != SWOverhangStrategy.SOFTCLIP ) {

            // we need to go through all the hassle below only if we do not do soft-clipping;
            // otherwise offset is never negative
            if ( offset < 0 ) {
                for (  ; j < (-offset) ; j++ ) {
                    bread.append((char) read[j]);
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
            for (; i < alignment.getAlignmentOffset() ; i++ ) {
                bref.append((char) ref[i]);
                bread.append(' ');
                match.append(' ');
            }
        }

        for ( final CigarElement e : cigar.getCigarElements() ) {
            switch (e.getOperator()) {
                case M :
                    for ( int z = 0 ; z < e.getLength() ; z++, i++, j++  ) {
                        bref.append((i< ref.length)?(char) ref[i]:' ');
                        bread.append((j < read.length)?(char) read[j]:' ');
                        match.append( ( i< ref.length && j < read.length ) ? (ref[i] == read[j] ? '.':'*' ) : ' ' );
                    }
                    break;
                case EQ:
                case X:
                    for ( int z = 0 ; z < e.getLength() ; z++, i++, j++  ) {
                        bref.append((i< ref.length)?(char) ref[i]:' ');
                        bread.append((j < read.length)?(char) read[j]:' ');
                        match.append( ( i< ref.length && j < read.length ) ? (e.getOperator()== CigarOperator.EQ ? '.':'*' ) : ' ' );
                    }
                    break;
                case I :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append('-');
                        bread.append((char) read[j]);
                        match.append('I');
                    }
                    break;
                case S :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append(' ');
                        bread.append((char) read[j]);
                        match.append('S');
                    }
                    break;
                case D:
                    for ( int z = 0 ; z < e.getLength(); z++ , i++ ) {
                        bref.append((char) ref[i]);
                        bread.append('-');
                        match.append('D');
                    }
                    break;
                default:
                    throw new GATKException("Unexpected Cigar element:" + e.getOperator());
            }
        }
        for (; i < ref.length; i++ ) bref.append((char) ref[i]);
        for (; j < read.length; j++ ) bread.append((char) read[j]);

        int pos = 0 ;
        final int maxlength = Math.max(match.length(), Math.max(bread.length(), bref.length()));
        while ( pos < maxlength ) {
            printCautiously(match, pos, 100);
            printCautiously(bread, pos, 100);
            printCautiously(bref, pos, 100);
            System.out.println();
            pos += 100;
        }
    }

    /** String builder's substring is extremely stupid: instead of trimming and/or returning an empty
     * string when one end/both ends of the interval are out of range, it crashes with an
     * exception. This utility function simply prints the substring if the interval is within the index range
     * or trims accordingly if it is not.
     */
    private static void printCautiously(final StringBuilder s, final int start, final int width) {
        if ( start >= s.length() ) {
            System.out.println();
            return;
        }
        final int end = Math.min(start + width, s.length());
        System.out.println(s.substring(start,end));
    }
}
