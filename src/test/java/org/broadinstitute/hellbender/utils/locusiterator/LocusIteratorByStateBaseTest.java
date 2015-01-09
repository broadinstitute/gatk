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

package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.BaseTest;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import java.util.*;
import java.util.stream.Collectors;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public class LocusIteratorByStateBaseTest extends BaseTest {
    protected static SAMFileHeader header;
    protected GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    private boolean isIndel(final CigarElement ce) {
        return ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I;
    }

    private boolean startsWithDeletion(final List<CigarElement> elements) {
        for ( final CigarElement element : elements ) {
            switch ( element.getOperator() ) {
                case M:
                case I:
                case EQ:
                case X:
                    return false;
                case D:
                    return true;
                default:
                    // keep looking
            }
        }

        return false;
    }

    private LIBSTest makePermutationTest(final List<CigarElement> elements) {
        CigarElement last = null;
        boolean hasMatch = false;

        // starts with D => bad
        if ( startsWithDeletion(elements) )
            return null;

        // ends with D => bad
        if ( elements.get(elements.size()-1).getOperator() == CigarOperator.D )
            return null;

        // make sure it's valid
        String cigar = "";
        for ( final CigarElement ce : elements ) {
            if ( ce.getOperator() == CigarOperator.N )
                return null; // TODO -- don't support N

            // abort on a bad cigar
            if ( last != null ) {
                if ( ce.getOperator() == last.getOperator() )
                    return null;
                if ( isIndel(ce) && isIndel(last) )
                    return null;
            }

            cigar += ce.getLength() + ce.getOperator().toString();
            last = ce;
            hasMatch = hasMatch || ce.getOperator() == CigarOperator.M;
        }

        if ( ! hasMatch && elements.size() == 1 &&
                ! (last.getOperator() == CigarOperator.I || last.getOperator() == CigarOperator.S))
            return null;

        return new LIBSTest(cigar);
    }

    @DataProvider(name = "LIBSTest")
    public Object[][] createLIBSTests(final List<Integer> cigarLengths, final List<Integer> combinations) {
        final List<Object[]> tests = new LinkedList<>();

        final List<CigarOperator> allOps = Arrays.asList(CigarOperator.values());

        final List<CigarElement> singleCigars = new LinkedList<>();
        for ( final int len : cigarLengths )
            singleCigars.addAll(allOps.stream().map(op -> new CigarElement(len, op)).collect(Collectors.toList()));

        for ( final int complexity : combinations ) {
            for ( final List<CigarElement> elements : Utils.makePermutations(singleCigars, complexity, true) ) {
                final LIBSTest test = makePermutationTest(elements);
                if ( test != null ) tests.add(new Object[]{test});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Work around inadequate tests that aren't worth fixing.
     *
     * Look at the CIGAR 2M2P2D2P2M.  Both M states border a deletion, separated by P (padding elements).  So
     * the right answer for deletions here is true for isBeforeDeletion() and isAfterDeletion() for the first
     * and second M.  But the LIBS_position doesn't say so.
     *
     * @param elements
     * @return
     */
    protected static boolean hasNeighboringPaddedOps(final List<CigarElement> elements, final int elementI) {
        return (elementI - 1 >= 0 && isPadding(elements.get(elementI-1))) ||
                (elementI + 1 < elements.size() && isPadding(elements.get(elementI+1)));
    }

    private static boolean isPadding(final CigarElement elt) {
        return elt.getOperator() == CigarOperator.P || elt.getOperator() == CigarOperator.H || elt.getOperator() == CigarOperator.S;
    }
}
