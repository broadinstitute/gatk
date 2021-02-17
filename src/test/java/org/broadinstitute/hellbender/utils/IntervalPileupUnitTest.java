package org.broadinstitute.hellbender.utils;

import com.google.common.collect.Iterators;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Unit tests for {@link IntervalPileup} implementations.
 */
public class IntervalPileupUnitTest extends GATKBaseTest {

    private static final String TEST_REFERENCE_FILE = b37_reference_20_21;
    private static final String TEST_ALIGNMENT_FILE = NA12878_20_21_WGS_bam;
    private static final SimpleInterval[] FIXED_TEST_INTERVALS = Arrays.stream(new String[] {"20:9,999,900-10,000,100", "21:9,999,900-10,000,100"})
            .map(SimpleInterval::new)
            .toArray(SimpleInterval[]::new);

    @DataProvider(name="sampledAndFixedIntervalData")
    public Iterator<Object[]> randomAndFixedIntervalData() {
        final Iterator<Object[]> fixed = Arrays.stream(fixedIntervalData()).iterator();
        final Iterator<Object[]> sampled = sampledData();

        return Iterators.concat(fixed, sampled);
    }

    @DataProvider(name = "sampledData")
    public Iterator<Object[]> sampledData() {
        final Path reference = IOUtils.getPath(TEST_REFERENCE_FILE);
        final int TELOMERE_SKIP = 20_000;
        final int SITES_PER_SEQ = 100;
        final ReferenceDataSource ref = new ReferenceFileSource(reference);
        final List<Object[]> result = new ArrayList<>(SITES_PER_SEQ * ref.getSequenceDictionary().getSequences().size());
        final int[] MARGINS = new int[] {5, 50, 100, 300}; // alternate with this.
        for (final SAMSequenceRecord sequenceRecord : ref.getSequenceDictionary().getSequences())  {
            final int spacing = (sequenceRecord.getSequenceLength() - 2 * TELOMERE_SKIP) / SITES_PER_SEQ;
            for (int i = 0; i < SITES_PER_SEQ; i++) {
                final int position = TELOMERE_SKIP + i * spacing;
                final int margin = MARGINS[i % MARGINS.length];
                final int start = Math.max(1, position - margin); // it should not go below 1 due to the telomere_skip but just being paranoic.
                final int end = Math.min(sequenceRecord.getSequenceLength(), position + margin);
                final SimpleInterval interval = new SimpleInterval(sequenceRecord.getSequenceName(), start, end);
                final ReadsDataSource readsDataSource = new ReadsPathDataSource(new GATKPath(TEST_ALIGNMENT_FILE));
                final List<GATKRead> reads = Utils.stream(readsDataSource.query(interval)).collect(Collectors.toList());
                final ReferenceBases refBases = new ReferenceBases(ref.queryAndPrefetch(interval).getBases(), interval);
                result.add(new Object[] { interval, reads, refBases});
            }
        }
        return result.iterator();
    }

    @Test(dataProvider = "sampledAndFixedIntervalData")
    public void testInterval(@SuppressWarnings("unused")
                             final SimpleInterval interval,
                             final List<GATKRead> reads,
                             final ReferenceBases referenceBases) {
        if (reads.size() > 2000) { // avoid test case that will take too long due to a large number of reads.
            Collections.shuffle(reads, new Random(reads.get(0).getName().hashCode()));
            reads.subList(2000, reads.size()).clear();
        }
        final TestIntervalPileup expected = new TestIntervalPileup(reads, referenceBases);
        final IntervalPileup actual = IntervalPileup.of(reads, referenceBases);
        Assert.assertEquals(actual.height(), expected.height());
        Assert.assertEquals(actual.width(), expected.width());
        Assert.assertEquals(actual.width(), referenceBases.getInterval().size());
        for (int i = 0; i < actual.height(); i++) {
            for (int j = 0; j < actual.width(); j++) {
                try {
                    Assert.assertEquals(actual.baseAt(i, j), expected.baseAt(i, j),
                            String.format("(%d, %d) (%d, %d)", i, j, (int) actual.baseAt(i, j), (int) expected.baseAt(i, j) ));

                    Assert.assertEquals(actual.qualAt(i, j), expected.qualAt(i, j),
                            String.format("(%d, %d)", i, j));
                    Assert.assertEquals(actual.hasInsertAt(i,j), expected.hasInsertAt(i,j),String.format("(%d, %d)", i, j));
                    assertEquals(actual.insertAt(i, j), expected.insertAt(i, j));
                    if (actual.hasInsertAt(i, j)) {
                        Assert.assertEquals(actual.insertBasesAt(i, j), expected.insertBasesAt(i, j), String.format("(%d, %d)", i, j));
                        Assert.assertEquals(actual.insertQualsAt(i, j), expected.insertQualsAt(i, j), String.format("(%d, %d)", i, j));
                    }
                } catch (final Error er) {
                    actual.insertBasesAt(i, j);
                    expected.insertBasesAt(i, j);
                    throw er;
                }
            }
        }
    }

    private void assertEquals(final IntervalPileup.Insert a, final IntervalPileup.Insert b) {
        if (a == b) return;
        Assert.assertEquals(a == null, b == null);
        Assert.assertEquals(a.length(), b.length());
        Assert.assertEquals(a.bases(), b.bases());
        Assert.assertEquals(a.quals(), b.quals());
        // extra checks.

        final byte[] aa = new byte[a.length() << 1];
        final byte[] bb = new byte[b.length()];
        a.copyBases(0, aa, 0, aa.length);
        b.copyBases(0, bb, 0, bb.length);
        Assert.assertTrue(Arrays.equals(Arrays.copyOf(aa, bb.length), bb));
        Assert.assertTrue(Arrays.equals(Arrays.copyOf(aa, bb.length), b.bases()));


        a.copyQuals(0, aa, 0, aa.length);
        b.copyQuals(0, bb, 0, bb.length);
        Assert.assertTrue(Arrays.equals(Arrays.copyOf(aa, bb.length), bb));
        Assert.assertTrue(Arrays.equals(Arrays.copyOf(aa, bb.length), b.quals()));
    }

    @DataProvider(name="fixedIntervalData")
    public Object[][] fixedIntervalData() {
        final Path reference = IOUtils.getPath(TEST_REFERENCE_FILE);
        final ReferenceDataSource ref = new ReferenceFileSource(reference);
        final ReadsDataSource readsDataSource = new ReadsPathDataSource(new GATKPath(TEST_ALIGNMENT_FILE));
        return Arrays.stream(FIXED_TEST_INTERVALS)
                .map(si -> new Object[] {si, Utils.stream(readsDataSource.query(si)).collect(Collectors.toList()),  new ReferenceBases(ref.queryAndPrefetch(si).getBases(), si)})
                .toArray(Object[][]::new);
    }

}

/**
 * Quicker and dirtier implementation of the interval pileup using a simpler and less efficient code.
 * <p>
 *     We will assume this implementation to be correct and we assert the run-time ones comparing the results
 *     using this implementation.
 * </p>
 * <p>
 * No code was borrowed from the main implementation so that mutually silencing bugs should be quite rare, at least
 * good for a prototype test.
 * </p>
 * <p>
 *     In order to make sure we will need to create controlled alignments
 * </p>
 */
class TestIntervalPileup implements IntervalPileup {

    private List<GATKRead> reads;

    private ReferenceBases referenceBases;

    TestIntervalPileup(final List<GATKRead> reads, final ReferenceBases ref) {
        this.reads = reads.stream()
                .filter(r -> !r.isUnmapped())
                .filter(r -> new SimpleInterval(r).overlaps(ref.getInterval()))
                .sorted(Comparator.comparingInt(GATKRead::getStart).thenComparing(GATKRead::getEnd).thenComparing(GATKRead::getName))
                .collect(Collectors.toList());
        this.referenceBases = ref;
    }

    @Override
    public List<GATKRead> reads() {
        return reads;
    }

    @Override
    public ReferenceBases reference() {
        return referenceBases;
    }

    @Override
    public int width() {
        return referenceBases.getInterval().size();
    }

    @Override
    public int height() {
        return reads.size();
    }

    @Override
    public byte baseAt(int row, int column) {
        final GATKRead read = reads.get(row);
        final int refPosition = column + referenceBases.getInterval().getStart();
        if (refPosition < read.getStart()) {
            return (byte) -1;
        } else if (refPosition > read.getEnd()) {
            return (byte) -1;
        } else {
            final int baseOffset =  matchReferencePosition(read.getCigar(), read.getStart(), refPosition);
            return baseOffset < 0 ? ((byte)'-') : read.getBase(baseOffset);
        }
    }

    @Override
    public boolean hasInsertAt(final int row, final int column) {
        final GATKRead read = reads.get(row);
        final Element element = element(read);
        return element.hasInsertAt(column);
    }

    @Override
    public byte[] insertBasesAt(final int row, final int column) {
        final GATKRead read = reads.get(row);
        final Element element = element(read);
        return element.insertBasesAt(column);
    }

    @Override
    public byte[] insertQualsAt(final int row, final int column) {
        final GATKRead read = reads.get(row);
        final Element element = element(read);
        return element.insertQualsAt(column);
    }


    @Override
    public byte qualAt(int row, int column) {
        final GATKRead read = reads.get(row);
        final int refPosition = column + referenceBases.getInterval().getStart();
        if (refPosition < read.getStart()) {
            return (byte) -1;
        } else if (refPosition > read.getEnd()) {
            return (byte) -1;
        } else {
            final int baseOffset =  matchReferencePosition(read.getCigar(), read.getStart(), refPosition);
            return baseOffset < 0 || ! read.hasBaseQualities()? ((byte)-1) : read.getBaseQuality(baseOffset);
        }
    }

    private int matchReferencePosition(final Cigar cigar, final int start, final int refPosition) {
        int i;
        int result = 0;
        for (i = 0; i < cigar.numCigarElements(); i++) {
            if (cigar.getCigarElement(i).getOperator().isAlignment()) {
                break;
            } else {
                result += cigar.getCigarElement(i).getOperator().consumesReadBases() ? cigar.getCigarElement(i).getLength() : 0;
            }
        }
        int refOffset = start;
        for (; i < cigar.numCigarElements(); i++) {
            final int len = cigar.getCigarElement(i).getLength();
            if (cigar.getCigarElement(i).getOperator().consumesReferenceBases()) {
                final int afterOffset = refOffset + len;
                if (afterOffset > refPosition) {
                    if (cigar.getCigarElement(i).getOperator().consumesReadBases()) {
                        return result + (refPosition - refOffset);
                    } else {
                        return -1;
                    }
                }
                refOffset = afterOffset;
            }
            if (cigar.getCigarElement(i).getOperator().consumesReadBases()) {
                result += len;
            }
        }
        Assert.fail("unreachable code");
        return -1;
    }

    @Override
    public GATKRead readAt(int row, int column) {
        final GATKRead read = reads.get(row);
        final int refPosition = referenceBases.getInterval().getStart() + column;
        if (refPosition < read.getStart()) {
            return null;
        } else if (refPosition > read.getEnd()) {
            return null;
        } else  {
            return read;
        }
    }

    @Override
    public List<GATKRead> readsAt(int row, int column) {
        return readAt(row, column) == null ? Collections.emptyList() : Collections.singletonList(readAt(row, column));
    }

    @Override
    public Element element(final GATKRead read) {
        final int index = reads.indexOf(read);
        if (index < 0) {
            return null;
        } else {
            return new TestElement(this, index, read);
        }
    }

class TestElement implements IntervalPileup.Element {

    private final int index;
    private final GATKRead read;
    private final TestIntervalPileup pileup;
    public TestElement(TestIntervalPileup pileup, int index, GATKRead read) {
        this.index = index;
        this.read = read;
        this.pileup = pileup;
    }

    @Override
    public GATKRead read() {
        return read;
    }

    @Override
    public int row() {
        return index;
    }

    @Override
    public int minColumn() {
        return Math.max(0, read.getStart() - referenceBases.getInterval().getStart());
    }

    @Override
    public int maxColumn() {
        return Math.max(referenceBases.getInterval().size(), read.getEnd() - referenceBases.getInterval().getStart() );
    }

    @Override
    public Insert insertAt(final int column) {
        final long search = findInsertStart(column);
        if (search == 0) {
            return null;
        } else {
            final byte[] bases = Arrays.copyOfRange(read.getBasesNoCopy(), (int) (search >> 32), ((int) (search >> 32)) + (int) search);
            final byte[] quals = Arrays.copyOfRange(read.getBaseQualitiesNoCopy(), (int) (search >> 32), ((int) (search >> 32)) + (int) search);
            return new TestInsert(column, bases, quals);
        }
    }

    @Override
    public List<Insert> inserts() {
        if (!read.getCigar().containsOperator(CigarOperator.I)) {
            return Collections.emptyList();
        } else {
            final List<Insert> result = new ArrayList<>(10);
            int i;
            final Cigar cigar = read.getCigar();
            int readOffset = 0;
            for (i = 0; i < cigar.numCigarElements(); i++) {
                if (cigar.getCigarElement(i).getOperator().isAlignment()) {
                    break;
                } else {
                    readOffset += cigar.getCigarElement(i).getOperator().consumesReadBases() ? cigar.getCigarElement(i).getLength() : 0;
                }
            }
            int refOffset = read.getStart();
            if (i > 0 && cigar.getCigarElement(i - 1).getOperator() == CigarOperator.I) {
                final int len = cigar.getCigarElement(i - 1).getLength();
                result.add(new TestInsert(refOffset, Arrays.copyOfRange(read.getBasesNoCopy(), readOffset, readOffset + len),
                        Arrays.copyOfRange(read.getBaseQualitiesNoCopy(), readOffset, readOffset + len)));
            }
            for (; i < cigar.numCigarElements(); i++) {
                if (cigar.getCigarElement(i).getOperator().consumesReferenceBases()) {
                    refOffset += cigar.getCigarElement(i).getLength();
                } else if (cigar.getCigarElement(i).getOperator() == CigarOperator.I) {
                    final int len = cigar.getCigarElement(i).getLength();
                    result.add(new TestInsert(refOffset, Arrays.copyOfRange(read.getBasesNoCopy(), readOffset, readOffset + len),
                            Arrays.copyOfRange(read.getBaseQualitiesNoCopy(), readOffset, readOffset + len)));
                }
                if (cigar.getCigarElement(i).getOperator().consumesReadBases()) {
                    readOffset += cigar.getCigarElement(i).getLength();
                }
            }
            return result;
        }
    }

    @Override
    public boolean hasInsertAt(int column) {
        return findInsertStart(column) != 0;
    }

    @Override
    public int insertSize(int column) {
        return (int) findInsertStart(column);
    }

    /**
     * Returns the offset within the read and the length that are typed as ints using
     * a long return. If there is no insertion (i.e. the length is 0) the the return
     * is simply zero.
     * If not zero, the lower 32 bits is the length of the insert, and the higher 32 bits the offset
     * within the read sequence.
     */
    private long findInsertStart(final int column) {
        final int where = column + referenceBases.getInterval().getStart();
        int i;
        final Cigar cigar = read.getCigar();
        int readOffset = 0;
        for (i = 0; i < cigar.numCigarElements(); i++) {
            if (cigar.getCigarElement(i).getOperator().isAlignment()) {
                break;
            } else {
                readOffset += cigar.getCigarElement(i).getOperator().consumesReadBases() ? cigar.getCigarElement(i).getLength() : 0;
            }
        }
        int refOffset = read.getStart();
        if (refOffset == where + 1 && i > 0 && cigar.getCigarElement(i - 1).getOperator() == CigarOperator.I) {
            final int len = cigar.getCigarElement(i - 1).getLength();
            return (((long)(readOffset - len)) << 32) + len;
        }
        for (; i < cigar.numCigarElements() && refOffset <= where + 1; i++) {
            if (cigar.getCigarElement(i).getOperator().consumesReferenceBases()) {
                refOffset += cigar.getCigarElement(i).getLength();
            } else if (cigar.getCigarElement(i).getOperator() == CigarOperator.I && refOffset == where + 1) {
                return (((long)(readOffset)) << 32) + cigar.getCigarElement(i).getLength();
            }
            if (cigar.getCigarElement(i).getOperator().consumesReadBases()) {
                readOffset += cigar.getCigarElement(i).getLength();
            }
        }
        return 0L;
    }

    @Override
    public int copyInsertBases(int column, byte[] dest, int offset, int maxLength) {
        final long search = findInsertStart(column);
        if (search == 0) {
            return 0;
        } else {
            final int length = (int) search;
            final int readOffset = (int) (search >> 32);
            return read.copyBases(readOffset, dest, offset, Math.min(maxLength, length));
        }
    }

    @Override
    public int copyInsertQuals(int column, byte[] dest, int offset, int maxLength) {
        final long search = findInsertStart(column);
        if (search == 0) {
            return 0;
        } else {
            final int length = (int) search;
            final int readOffset = (int) (search >> 32);
            return read.copyBaseQualities(readOffset, dest, offset, Math.min(maxLength, length));
        }
    }

    @Override
    public byte[] insertQualsAt(int column) {
        final int insertLength = insertSize(column);
        final byte[] result = new byte[insertLength];
        final int copied = copyInsertQuals(column, result, 0, insertLength);
        return copied == insertLength ? result : Arrays.copyOf(result, copied);
    }

    @Override
    public byte[] insertBasesAt(int column) {
        final int insertLength = insertSize(column);
        final byte[] result = new byte[insertLength];
        final int copied = copyInsertBases(column, result, 0, insertLength);
        return copied == insertLength ? result : Arrays.copyOf(result, copied);
    }

    @Override
    public byte baseAt(int column) {
        return pileup.baseAt(index, column);
    }

    @Override
    public byte qualAt(int column) {
        return pileup.qualAt(index, column);
    }
    }

    @Override
    public Insert insertAt(final int row, final int column) {
        final GATKRead read = reads.get(row);
        final Element element = new TestElement(this, row, read);
        return element.insertAt(column);
    }
}

class TestInsert implements IntervalPileup.Insert {

    private final byte[] bases;
    private final byte[] quals;
    private final int column;

    TestInsert(final int column, final byte[] bases, final byte[] quals) {
        this.column = column;
        this.bases = bases;
        this.quals = quals;
    }

    public int column() {
        return column;
    }

    @Override
    public int length() {
        return bases.length;
    }

    @Override
    public byte[] bases() {
        return bases;
    }

    @Override
    public byte[] quals() {
        return quals;
    }

    @Override
    public int copyBases(int offset, byte[] dest, int destOffset, int maxLength) {
        final int length = Math.min(bases.length - offset, maxLength);
        System.arraycopy(bases, offset, dest, destOffset, length);
        return length;
    }

    @Override
    public int copyQuals(int offset, byte[] dest, int destOffset, int maxLength) {
        final int length = Math.min(quals.length - offset, maxLength);
        System.arraycopy(quals, offset, dest, destOffset, length);
        return length;
    }

    public int hashCode() {
        return Arrays.hashCode(bases);
    }

    @Override
    public boolean equals(final Object other) {
        if (other instanceof IntervalPileup.Insert) {
            final IntervalPileup.Insert insert = (IntervalPileup.Insert) other;
            return Arrays.equals(bases, insert.bases()) && Arrays.equals(quals, insert.quals());
        } else {
            return false;
        }
    }

    @Override
    public String toString() {
            return new String(bases()) + "/" + new String(quals());
    }
}