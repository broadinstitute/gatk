package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;

import java.util.*;
import java.util.stream.Collectors;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public abstract class LocusIteratorByStateBaseTest extends GATKBaseTest {
    protected static SAMFileHeader header;
    protected GenomeLocParser genomeLocParser;

    /**
     * For testing only.  Assumes that the incoming SAMRecords have no read groups, so creates a dummy sample list
     * for the system.
     */
    public static List<String> sampleListForSAMWithoutReadGroups() {
        List<String> samples = new ArrayList<>();
        samples.add(null);
        return samples;
    }

    @BeforeClass
    public void beforeClass() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    protected LocusIteratorByState makeLIBS(final List<GATKRead> reads, final SAMFileHeader header) {
        return makeLIBS(reads, null, false, header);
    }

    protected LocusIteratorByState makeLIBSwithNs(final List<GATKRead> reads, final SAMFileHeader header) {
        return makeLIBS(reads, null, true, false, header);
    }

    protected LocusIteratorByState makeLIBS(final List<GATKRead> reads,
        final DownsamplingMethod downsamplingMethod,
        final boolean keepUniqueReadList,
        final SAMFileHeader header) {
        return makeLIBS(reads, downsamplingMethod, false, keepUniqueReadList, header);
    }

    protected LocusIteratorByState makeLIBS(final List<GATKRead> reads,
                                            final DownsamplingMethod downsamplingMethod,
                                            final boolean includeNs,
                                            final boolean keepUniqueReadList,
                                            final SAMFileHeader header) {
        reads.sort(new ReadCoordinateComparator(header));
        return new LocusIteratorByState(
                new FakeCloseableIterator<>(reads.iterator()),
                downsamplingMethod,
                keepUniqueReadList,
                sampleListForSAMWithoutReadGroups(),
                header,
                true,
                includeNs
        );
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


    static final class FakeCloseableIterator<T> implements CloseableIterator<T> {
        Iterator<T> iterator;

        public FakeCloseableIterator(Iterator<T> it) {
            iterator = it;
        }

        @Override
        public void close() {}

        @Override
        public boolean hasNext() {
            return iterator.hasNext();
        }

        @Override
        public T next() {
            return iterator.next();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Don't remove!");
        }
    }
}
