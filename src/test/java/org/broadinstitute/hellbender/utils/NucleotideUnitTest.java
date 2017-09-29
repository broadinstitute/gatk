package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for {@link Nucleotide}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class NucleotideUnitTest {

    private static final Random random = new Random(13);
    private static final RandomDNA randomDNA = new RandomDNA(13);
    private static final int MIN_RANDOM_SEQ_LENGTH = 10;
    private static final int MAX_RANDOM_SEQ_LENGTH = 100;
    private static final int NUMBER_OF_RANDOM_SEQUENCES = 10;

    @Test
    public void testToBase() {
        Assert.assertEquals(Nucleotide.A.toBase(), (byte)'A');
        Assert.assertEquals(Nucleotide.C.toBase(), (byte)'C');
        Assert.assertEquals(Nucleotide.G.toBase(), (byte)'G');
        Assert.assertEquals(Nucleotide.N.toBase(), (byte)'N');
        Assert.assertEquals(Nucleotide.T.toBase(), (byte)'T');
        Assert.assertEquals(Nucleotide.X.toBase(), (byte)'X');
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testToBaseOnInvalid() {
        Nucleotide.INVALID.toBase();
    }

    @Test
    public void testValueOfBase() {
        for (byte i = 0; i >= 0; i++) {
            final Nucleotide expected;
            switch (i) {
                case 'a':
                case 'A': expected = Nucleotide.A; break;
                case 'c':
                case 'C': expected = Nucleotide.C; break;
                case 'g':
                case 'G': expected = Nucleotide.G; break;
                case 't':
                case 'T':
                case 'u':
                case 'U': expected = Nucleotide.T; break;
                case 'n':
                case 'N': expected = Nucleotide.N; break;
                case 'x':
                case 'X': expected = Nucleotide.X; break;
                default : expected = Nucleotide.INVALID;
            }
            Assert.assertSame(Nucleotide.valueOf(i), expected, "Failed with base " + i + " returning nucleotide " + Nucleotide.valueOf(i));
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testValueOfNegativeBase() {
        Nucleotide.valueOf((byte) -10);
    }

    @Test
    public void testNucleotideCounterInit() {
        final Nucleotide.Counter counter = new Nucleotide.Counter();
        for (final Nucleotide n : Nucleotide.values()) {
            Assert.assertEquals(counter.get(n), 0);
        }
    }

    @Test(dependsOnMethods = "testValueOfBase", dataProvider = "testSequences")
    public void testAddingOneByOne(final byte[] bases) {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        final Map<Nucleotide, Integer> shadow = new HashMap<>(Nucleotide.values().length);
        for (final byte base : bases) {
            subject.add(base);
            final Nucleotide nuc = Nucleotide.valueOf(base);
            shadow.put(nuc, shadow.getOrDefault(nuc, 0) + 1);
            for (final Nucleotide n : Nucleotide.values()) {
                Assert.assertEquals(subject.get(n), (long) shadow.getOrDefault(n, 0));
            }
        }
    }

    @Test(dependsOnMethods = "testValueOfBase", dataProvider = "testSequences")
    public void testAddingAllAtOnce(final byte[] bases) {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        final Map<Nucleotide, Integer> shadow = new HashMap<>(Nucleotide.values().length);
        for (final byte base : bases) {
            final Nucleotide nuc = Nucleotide.valueOf(base);
            shadow.put(nuc, shadow.getOrDefault(nuc, 0) + 1);
        }
        subject.addAll(bases);
        for (final Nucleotide n : Nucleotide.values()) {
            Assert.assertEquals(subject.get(n), (long) shadow.getOrDefault(n, 0));
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddingAllAtOnceOnANullArray() {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        subject.addAll((byte[])null);
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddingAllAtOnceWithNegativeBases() {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        subject.addAll(new byte[] { 'a', 'A', -10, 'C' } );
    }



    @Test(dependsOnMethods = "testValueOfBase", dataProvider = "testSequences")
    public void testClear(final byte[] bases) {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        final Map<Nucleotide, Integer> shadow = new HashMap<>(Nucleotide.values().length);
        for (final byte base : bases) {
            final Nucleotide nuc = Nucleotide.valueOf(base);
            shadow.put(nuc, shadow.getOrDefault(nuc, 0) + 1);
        }
        subject.addAll(bases);
        for (final Nucleotide n : Nucleotide.values()) {
            Assert.assertEquals(subject.get(n), (long) shadow.getOrDefault(n, 0));
        }
        subject.clear();
        for (final Nucleotide n : Nucleotide.values()) {
            Assert.assertEquals(subject.get(n), 0);
        }
    }

    @DataProvider(name = "testSequences")
    public Object[][] testSequences() {
        final List<Object[]> result = new ArrayList<>();
        // We add non random trivial sequences:
        result.add(new Object[] { new byte[0] });
        for (final Nucleotide nuc : Nucleotide.values()) {
            if (nuc == Nucleotide.INVALID) {
                continue;
            }
            result.add( new Object[] { new byte[] { nuc.toBase() } });
            result.add( new Object[] { Utils.repeatBytes( nuc.toBase(), MIN_RANDOM_SEQ_LENGTH) });
        }
        for (int i = 0; i < NUMBER_OF_RANDOM_SEQUENCES; i++) {
            final int length = random.nextInt(MAX_RANDOM_SEQ_LENGTH - MIN_RANDOM_SEQ_LENGTH + 1) + MIN_RANDOM_SEQ_LENGTH;
            final byte[] base = randomDNA.nextBases(length);
            result.add(new Object[] { base });
        }
        return result.toArray(new Object[result.size()][]);
    }
}
