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
        Assert.assertEquals(Nucleotide.A.encodeAsByte(), (byte)'A');
        Assert.assertEquals(Nucleotide.C.encodeAsByte(), (byte)'C');
        Assert.assertEquals(Nucleotide.G.encodeAsByte(), (byte)'G');
        Assert.assertEquals(Nucleotide.N.encodeAsByte(), (byte)'N');
        Assert.assertEquals(Nucleotide.T.encodeAsByte(), (byte)'T');
        Assert.assertEquals(Nucleotide.X.encodeAsByte(), (byte)'X');
    }

    @Test
    public void testIsConcrete() {
        for (final Nucleotide nuc : Nucleotide.values()) {
            switch (nuc) {
                case A:
                case C:
                case T:
                case G:
                    Assert.assertTrue(nuc.isConcrete());
                    break;
                default:
                    Assert.assertFalse(nuc.isConcrete());
            }
        }
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testToBaseOnInvalid() {
        Nucleotide.INVALID.encodeAsByte();
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
                case 'r':
                case 'R': expected = Nucleotide.R; break;
                case 'b':
                case 'B': expected = Nucleotide.B; break;
                case 'v':
                case 'V': expected = Nucleotide.V; break;
                case 'y':
                case 'Y': expected = Nucleotide.Y; break;
                case 's':
                case 'S': expected = Nucleotide.S; break;
                case 'w':
                case 'W': expected = Nucleotide.W; break;
                case 'k':
                case 'K': expected = Nucleotide.K; break;
                case 'm':
                case 'M': expected = Nucleotide.M; break;
                case 'd':
                case 'D': expected = Nucleotide.D; break;
                case 'h':
                case 'H': expected = Nucleotide.H; break;
                default :
                    expected = Nucleotide.INVALID;
            }
            Assert.assertSame(Nucleotide.decode(i), expected, "Failed with base " + i + " returning nucleotide " + Nucleotide.decode(i));
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testValueOfNegativeBase() {
        Nucleotide.decode((byte) -10);
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
            final Nucleotide nuc = Nucleotide.decode(base);
            shadow.put(nuc, shadow.getOrDefault(nuc, 0) + 1);
            for (final Nucleotide n : Nucleotide.values()) {
                Assert.assertEquals(subject.get(n), (long) shadow.getOrDefault(n, 0));
            }
        }
        Assert.assertEquals(subject.sum(), shadow.values().stream().mapToLong(l -> l).sum());
    }

    @Test(dependsOnMethods = "testValueOfBase", dataProvider = "testSequences")
    public void testAddingAllAtOnce(final byte[] bases) {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        final Map<Nucleotide, Integer> shadow = new HashMap<>(Nucleotide.values().length);
        for (final byte base : bases) {
            final Nucleotide nuc = Nucleotide.decode(base);
            shadow.put(nuc, shadow.getOrDefault(nuc, 0) + 1);
        }
        subject.addAll(bases);
        for (final Nucleotide n : Nucleotide.values()) {
            Assert.assertEquals(subject.get(n), (long) shadow.getOrDefault(n, 0));
        }
        Assert.assertEquals(subject.sum(), shadow.values().stream().mapToLong(l -> l).sum());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddingAllAtOnceOnANullArray() {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        subject.addAll(null);
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
            final Nucleotide nuc = Nucleotide.decode(base);
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
        Assert.assertEquals(subject.sum(), 0);
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
            result.add( new Object[] { new byte[] { nuc.encodeAsByte() } });
            result.add( new Object[] { Utils.repeatBytes( nuc.encodeAsByte(), MIN_RANDOM_SEQ_LENGTH) });
        }
        for (int i = 0; i < NUMBER_OF_RANDOM_SEQUENCES; i++) {
            final int length = random.nextInt(MAX_RANDOM_SEQ_LENGTH - MIN_RANDOM_SEQ_LENGTH + 1) + MIN_RANDOM_SEQ_LENGTH;
            final byte[] base = randomDNA.nextBases(length);
            result.add(new Object[] { base });
        }
        return result.toArray(new Object[result.size()][]);
    }
}
