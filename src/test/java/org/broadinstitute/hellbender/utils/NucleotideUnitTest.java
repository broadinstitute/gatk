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
    private static final RandomDNA randomDNA = new RandomDNA(random);
    private static final int MIN_RANDOM_SEQ_LENGTH = 10;
    private static final int MAX_RANDOM_SEQ_LENGTH = 100;
    private static final int NUMBER_OF_RANDOM_SEQUENCES = 10;

    @Test(dataProvider = "values")
    public void testEncodeAsByte(final Nucleotide nuc) {
        // Will always use the first letter of the constant as the one byte encoding.
        final char firstLetter = nuc.name().charAt(0);
        final byte expectedLowerEncoding = (byte) Character.toLowerCase(firstLetter);
        final byte expectedUpperEncoding = (byte) Character.toUpperCase(firstLetter);
        Assert.assertEquals(nuc.encodeAsByte(), expectedUpperEncoding); // by default is upper case.
        Assert.assertEquals(nuc.encodeAsByte(true), expectedUpperEncoding);
        Assert.assertEquals(nuc.encodeAsByte(false), expectedLowerEncoding);
    }

    @Test(dataProvider = "values")
    public void testEncodeAsChar(final Nucleotide nuc) {
        // Will always use the first letter of the constant as the one byte encoding.
        final char firstLetter = nuc.name().charAt(0);
        final char expectedLowerEncoding = Character.toLowerCase(firstLetter);
        final char expectedUpperEncoding = Character.toUpperCase(firstLetter);
        Assert.assertEquals(nuc.encodeAsChar(), expectedUpperEncoding); // by default is upper case.
        Assert.assertEquals(nuc.encodeAsChar(true), expectedUpperEncoding);
        Assert.assertEquals(nuc.encodeAsChar(false), expectedLowerEncoding);
    }

    @Test(dataProvider = "values")
    public void testEncodeAsString(final Nucleotide nuc) {
        Assert.assertEquals(nuc.encodeAsString(), "" + (char) nuc.encodeAsByte());
        Assert.assertEquals(nuc.encodeAsString(), nuc.encodeAsString(true));
        Assert.assertEquals(nuc.encodeAsString(false), nuc.encodeAsString(true).toLowerCase());
        Assert.assertEquals(nuc.encodeAsString(true), "" + (char) nuc.encodeAsByte(true));
        Assert.assertEquals(nuc.encodeAsString(false), ""  + (char) nuc.encodeAsByte(false));
    }

    /**
     * Checks the assumption that each nuc canonical name is a single upper-case letter.
     * If this test fails that would  indicates that you are modifying {@link Nucleotide} in a way
     * that may result in unexpected errors down-the-road as some code that uses that class
     * may make such an assumption.
     *
     * @param nuc the value to test.
     */
    @Test(dataProvider = "values")
    public void testSingleUpperLetterNames(final Nucleotide nuc) {
        Assert.assertEquals(nuc.name().length(), 1);
        Assert.assertEquals(nuc.name(), nuc.name().toUpperCase());
        Assert.assertEquals(nuc.encodeAsString().length(), 1);
        Assert.assertEquals(nuc.encodeAsString(false).length(), 1);
    }

    @Test
    public void testStandardNucleotidesList() {
        Assert.assertEquals(Nucleotide.STANDARD_BASES.size(), 4);
        for (final Nucleotide value : Nucleotide.values()) {
            switch (value) {
                case A:
                case C:
                case G:
                case T:
                    Assert.assertTrue(Nucleotide.STANDARD_BASES.contains(value));
                    break;
                default:
                    Assert.assertFalse(Nucleotide.STANDARD_BASES.contains(value));
            }
        }
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testStandardNucleotideListIsUnmodifiable1() {
        Nucleotide.STANDARD_BASES.clear();
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testStandardNucleotideListIsUnmodifiable2() {
        Nucleotide.STANDARD_BASES.set(0, Nucleotide.N);
    }

    @Test(dataProvider = "values")
    public void testIsStandard(final Nucleotide nuc) {
        switch (nuc) {
            case A:
            case C:
            case T:
            case G:
                Assert.assertTrue(nuc.isStandard());
                break;
            default:
                Assert.assertFalse(nuc.isStandard());
        }
    }

    @Test(dataProvider = "values")
    public void testIsAmbiguous(final Nucleotide nuc) {
        switch (nuc) {
            case X:
            case A:
            case C:
            case T:
            case G:
                Assert.assertFalse(nuc.isAmbiguous());
                break;
            default:
                Assert.assertTrue(nuc.isAmbiguous());
        }
    }

    @Test(dataProvider = "values")
    public void testIsValid(final Nucleotide nuc) {
        switch (nuc) {
            case X:
                Assert.assertFalse(nuc.isValid());
                break;
            default:
                Assert.assertTrue(nuc.isValid());
        }
    }

    @Test
    public void testDecode() {
        for (byte i = 0; i >= 0; i++) {
            final Nucleotide expected;
            switch (i) {
                case 'a':
                case 'A':
                    expected = Nucleotide.A;
                    break;
                case 'c':
                case 'C':
                    expected = Nucleotide.C;
                    break;
                case 'g':
                case 'G':
                    expected = Nucleotide.G;
                    break;
                case 't':
                case 'T':
                case 'u':
                case 'U':
                    expected = Nucleotide.T;
                    break;
                case 'n':
                case 'N':
                    expected = Nucleotide.N;
                    break;
                case 'x':
                case 'X':
                    expected = Nucleotide.X;
                    break;
                case 'r':
                case 'R':
                    expected = Nucleotide.R;
                    break;
                case 'b':
                case 'B':
                    expected = Nucleotide.B;
                    break;
                case 'v':
                case 'V':
                    expected = Nucleotide.V;
                    break;
                case 'y':
                case 'Y':
                    expected = Nucleotide.Y;
                    break;
                case 's':
                case 'S':
                    expected = Nucleotide.S;
                    break;
                case 'w':
                case 'W':
                    expected = Nucleotide.W;
                    break;
                case 'k':
                case 'K':
                    expected = Nucleotide.K;
                    break;
                case 'm':
                case 'M':
                    expected = Nucleotide.M;
                    break;
                case 'd':
                case 'D':
                    expected = Nucleotide.D;
                    break;
                case 'h':
                case 'H':
                    expected = Nucleotide.H;
                    break;
                default:
                    expected = Nucleotide.X;
            }
            Assert.assertSame(Nucleotide.decode(i), expected, "Failed with base " + i + " returning nucleotide " + Nucleotide.decode(i));
            Assert.assertSame(Nucleotide.decode((char)i), expected, "Failed with base " + i + " returning nucleotide " + Nucleotide.decode((char)i));
            final StringBuilder builder = new StringBuilder(1);
            builder.append((char)i);
            Assert.assertSame(Nucleotide.decode(builder.toString()), expected);
        }
    }

    @Test
    public void testDecodeCharOverMax() {
        final char baseCh = 'A';
        final char alteredCh = (char) (baseCh | 0x0100);
        Assert.assertSame(Nucleotide.decode(baseCh), Nucleotide.A);
        Assert.assertSame(Nucleotide.decode(alteredCh), Nucleotide.X);
    }

    @Test(dataProvider = "values")
    public void testIncludes(final Nucleotide nuc) {
        if (nuc.isStandard()) {
            for (final Nucleotide other : Nucleotide.values()) {
                if (other.isStandard()) {
                    Assert.assertEquals(nuc.includes(other), nuc == other);
                    Assert.assertEquals(nuc.includes(other.encodeAsByte()), nuc == other);
                    Assert.assertEquals(nuc.includes(other.encodeAsByte(false)), nuc == other);
                } else {
                    Assert.assertFalse(nuc.includes(other));
                    Assert.assertFalse(nuc.includes(other.encodeAsByte()));
                    Assert.assertFalse(nuc.includes(other.encodeAsByte(false)));
                }
            }
        } else if (nuc.isAmbiguous()) {
            for (final Nucleotide other : Nucleotide.values()) {
                final boolean thisA = nuc.includes(Nucleotide.A);
                final boolean thisC = nuc.includes(Nucleotide.C);
                final boolean thisG = nuc.includes(Nucleotide.G);
                final boolean thisT = nuc.includes(Nucleotide.T);
                final boolean otherA = other.includes(Nucleotide.A);
                final boolean otherC = other.includes(Nucleotide.C);
                final boolean otherG = other.includes(Nucleotide.G);
                final boolean otherT = other.includes(Nucleotide.T);
                final boolean includes = other.isValid() && (thisA == otherA || thisA)
                        && (thisC == otherC || thisC)
                        && (thisG == otherG || thisG)
                        && (thisT == otherT || thisT);
                Assert.assertEquals(nuc.includes(other), includes, "" + nuc + " " + other);
                Assert.assertEquals(nuc.includes(other.encodeAsByte()), includes);
                Assert.assertEquals(nuc.includes(other.encodeAsByte(false)), includes);
            }
        } else { // invalid
            for (final Nucleotide other : Nucleotide.values()) {
                Assert.assertFalse(nuc.includes(other));
                Assert.assertFalse(nuc.includes(other.encodeAsByte()));
                Assert.assertFalse(nuc.includes(other.encodeAsByte(false)));
            }
        }
    }

    @Test(dataProvider = "values")
    public void testIntersects(final Nucleotide nuc) {
        final boolean thisA = nuc.includes(Nucleotide.A);
        final boolean thisC = nuc.includes(Nucleotide.C);
        final boolean thisG = nuc.includes(Nucleotide.G);
        final boolean thisT = nuc.includes(Nucleotide.T);
        for (final Nucleotide other : Nucleotide.values()) {
            final boolean otherA = other.includes(Nucleotide.A);
            final boolean otherC = other.includes(Nucleotide.C);
            final boolean otherG = other.includes(Nucleotide.G);
            final boolean otherT = other.includes(Nucleotide.T);
            final Nucleotide intersect = nuc.intersect(other);
            Assert.assertNotNull(intersect);
            Assert.assertEquals(intersect.includes(Nucleotide.A), thisA && otherA);
            Assert.assertEquals(intersect.includes(Nucleotide.C), thisC && otherC);
            Assert.assertEquals(intersect.includes(Nucleotide.G), thisG && otherG);
            Assert.assertEquals(intersect.includes(Nucleotide.T), thisT && otherT);
        }
    }

    @Test(dataProvider = "values")
    public void testComplement(final Nucleotide nuc) {
        final boolean thisA = nuc.includes(Nucleotide.A);
        final boolean thisC = nuc.includes(Nucleotide.C);
        final boolean thisG = nuc.includes(Nucleotide.G);
        final boolean thisT = nuc.includes(Nucleotide.T);
        final Nucleotide complement = nuc.complement();
        final boolean compA = complement.includes(Nucleotide.A);
        final boolean compC = complement.includes(Nucleotide.C);
        final boolean compG = complement.includes(Nucleotide.G);
        final boolean compT = complement.includes(Nucleotide.T);
        final String errorMessage = "Failure with " + nuc + " result in complement " + complement;
        Assert.assertEquals(compA, thisT, errorMessage);
        Assert.assertEquals(compT, thisA, errorMessage);
        Assert.assertEquals(compC, thisG, errorMessage);
        Assert.assertEquals(compG, thisC, errorMessage);
    }

    @Test(dataProvider = "values")
    public void testTransition(final Nucleotide nuc) {
        final boolean thisA = nuc.includes(Nucleotide.A);
        final boolean thisC = nuc.includes(Nucleotide.C);
        final boolean thisG = nuc.includes(Nucleotide.G);
        final boolean thisT = nuc.includes(Nucleotide.T);
        final Nucleotide trans = nuc.transition();
        final boolean tranA = trans.includes(Nucleotide.A);
        final boolean tranC = trans.includes(Nucleotide.C);
        final boolean tranG = trans.includes(Nucleotide.G);
        final boolean tranT = trans.includes(Nucleotide.T);
        final String errorMessage = "Failure with " + nuc + " result in transition " + trans;
        Assert.assertEquals(tranA, thisG, errorMessage);
        Assert.assertEquals(tranG, thisA, errorMessage);
        Assert.assertEquals(tranC, thisT, errorMessage);
        Assert.assertEquals(tranT, thisC, errorMessage);
    }

    @Test(dataProvider = "values")
    public void testTransversion(final Nucleotide nuc) {
        final boolean thisA = nuc.includes(Nucleotide.A);
        final boolean thisC = nuc.includes(Nucleotide.C);
        final boolean thisG = nuc.includes(Nucleotide.G);
        final boolean thisT = nuc.includes(Nucleotide.T);
        final Nucleotide trans = nuc.transversion();
        final boolean tranA = trans.includes(Nucleotide.A);
        final boolean tranC = trans.includes(Nucleotide.C);
        final boolean tranG = trans.includes(Nucleotide.G);
        final boolean tranT = trans.includes(Nucleotide.T);
        final String errorMessage = "Failure with " + nuc + " result in transversion " + trans;
        Assert.assertEquals(tranA, thisC || thisT, errorMessage);
        Assert.assertEquals(tranG, thisC || thisT, errorMessage);
        Assert.assertEquals(tranC, thisA || thisG, errorMessage);
        Assert.assertEquals(tranT, thisA || thisG, errorMessage);
        final Nucleotide transStrong = nuc.transversion(true);
        final Nucleotide transWeak = nuc.transversion(false);
        Assert.assertTrue(trans.includes(transStrong) || trans == Nucleotide.X || transStrong == Nucleotide.X);
        Assert.assertTrue(trans.includes(transWeak) || trans == Nucleotide.X || transStrong == Nucleotide.X);
        Assert.assertSame(transStrong.intersect(transWeak), Nucleotide.X);
        Assert.assertEquals(transStrong.includes(Nucleotide.C), tranC);
        Assert.assertEquals(transStrong.includes(Nucleotide.G), tranG);
        Assert.assertEquals(transWeak.includes(Nucleotide.A), tranA);
        Assert.assertEquals(transWeak.includes(Nucleotide.T), tranT);
    }

    @Test(dataProvider = "values")
    public void testSame(final Nucleotide nuc) {
        for (final Nucleotide other : Nucleotide.values()) {
            final boolean reallyTheSame = nuc != Nucleotide.INVALID && nuc == other;
            Assert.assertEquals(nuc.same(other), reallyTheSame);
            Assert.assertEquals(Nucleotide.same(nuc.encodeAsByte(), other.encodeAsByte()), reallyTheSame);
            Assert.assertEquals(Nucleotide.same(nuc.encodeAsByte(false), other.encodeAsByte()), reallyTheSame);
            Assert.assertEquals(Nucleotide.same(nuc.encodeAsByte(), other.encodeAsByte(false)), reallyTheSame);
            Assert.assertEquals(Nucleotide.same(nuc.encodeAsByte(false), other.encodeAsByte(false)), reallyTheSame);
        }
    }

    @Test
    public void testValueOfNegativeBase() {
        Assert.assertSame(Nucleotide.decode((byte) -10), Nucleotide.X);
    }

    @Test
    public void testNucleotideCounterInit() {
        final Nucleotide.Counter counter = new Nucleotide.Counter();
        for (final Nucleotide n : Nucleotide.values()) {
            Assert.assertEquals(counter.get(n), 0);
        }
    }

    @Test(dependsOnMethods = "testDecode", dataProvider = "testSequences")
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

    @Test(dependsOnMethods = {"testDecode", "testDecodeCharOverMax"}, dataProvider = "testSequences")
    public void testAddingCharsOneByOne(final byte[] byteBases) {
        final char[] bases = new char[byteBases.length];
        for (int i = 0; i < bases.length; i++) {
            bases[i] = (char) byteBases[i];
            if ((random.nextDouble() < 0.1)) {
                // 10% of the time we will add random fuzz to the upper-byte of the char to render it
                // invalid.
                bases[i] |= (random.nextInt(127) + 1) << 8;
            }
        }
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        final Map<Nucleotide, Integer> shadow = new HashMap<>(Nucleotide.values().length);
        for (final char base : bases) {
            subject.add(base);
            final Nucleotide nuc = Nucleotide.decode(base);
            shadow.put(nuc, shadow.getOrDefault(nuc, 0) + 1);
            for (final Nucleotide n : Nucleotide.values()) {
                Assert.assertEquals(subject.get(n), (long) shadow.getOrDefault(n, 0));
            }
        }
        Assert.assertEquals(subject.sum(), shadow.values().stream().mapToLong(l -> l).sum());
    }

    @Test(dependsOnMethods = "testDecode", dataProvider = "testSequences")
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

    @Test(dependsOnMethods = {"testDecode", "testDecodeCharOverMax"}, dataProvider = "testSequences")
    public void testAddingAllCharsAtOnce(final byte[] byteBases) {
        final char[] bases = new char[byteBases.length];
        for (int i = 0; i < bases.length; i++) {
            bases[i] = (char) byteBases[i];
            if ((random.nextDouble() < 0.1)) {
                // 10% of the time we will add random fuzz to the upper-byte of the char to render it
                // invalid.
                bases[i] |= (random.nextInt(127) + 1) << 8;
            }
        }
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        final Map<Nucleotide, Integer> shadow = new HashMap<>(Nucleotide.values().length);
        for (final char base : bases) {
            final Nucleotide nuc = Nucleotide.decode(base);
            shadow.put(nuc, shadow.getOrDefault(nuc, 0) + 1);
        }
        subject.addAll(bases);
        for (final Nucleotide n : Nucleotide.values()) {
            Assert.assertEquals(subject.get(n), (long) shadow.getOrDefault(n, 0));
        }
        Assert.assertEquals(subject.sum(), shadow.values().stream().mapToLong(l -> l).sum());
    }

    @Test(dependsOnMethods = {"testDecode", "testDecodeCharOverMax"}, dataProvider = "testSequences")
    public void testAddingStringAtOnce(final byte[] byteBases) {
        final char[] bases = new char[byteBases.length];
        for (int i = 0; i < bases.length; i++) {
            bases[i] = (char) byteBases[i];
            if ((random.nextDouble() < 0.1)) {
                // 10% of the time we will add random fuzz to the upper-byte of the char to render it
                // invalid.
                bases[i] |= (random.nextInt(127) + 1) << 8;
            }
        }
        final String string = new String(bases);
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        final Map<Nucleotide, Integer> shadow = new HashMap<>(Nucleotide.values().length);
        for (final char base : bases) {
            final Nucleotide nuc = Nucleotide.decode(base);
            shadow.put(nuc, shadow.getOrDefault(nuc, 0) + 1);
        }
        subject.addAll(string);
        for (final Nucleotide n : Nucleotide.values()) {
            Assert.assertEquals(subject.get(n), (long) shadow.getOrDefault(n, 0));
        }
        Assert.assertEquals(subject.sum(), shadow.values().stream().mapToLong(l -> l).sum());
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddingAllAtOnceOnANullArray() {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        subject.addAll((byte[]) null);
    }

    @Test
    public void testAddingAllAtOnceWithNegativeBases() {
        final Nucleotide.Counter subject = new Nucleotide.Counter();
        subject.addAll(new byte[]{'a', 'A', -10, 'C'});
        Assert.assertEquals(subject.get(Nucleotide.INVALID), 1);
        Assert.assertEquals(subject.get(Nucleotide.A), 2);
        Assert.assertEquals(subject.get(Nucleotide.T), 0);
    }

    @Test(dependsOnMethods = "testDecode", dataProvider = "testSequences")
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

    @Test
    public void testUracilSameAsThymine() {
        Assert.assertSame(Nucleotide.U, Nucleotide.T);
    }

    @Test
    public void testLongFormNames() {
        Assert.assertSame(Nucleotide.ADENINE, Nucleotide.A);
        Assert.assertSame(Nucleotide.THYMINE, Nucleotide.T);
        Assert.assertSame(Nucleotide.GUANINE, Nucleotide.G);
        Assert.assertSame(Nucleotide.CYTOSINE, Nucleotide.C);
        Assert.assertSame(Nucleotide.URACIL, Nucleotide.U);
        Assert.assertSame(Nucleotide.ANY, Nucleotide.N);
        Assert.assertSame(Nucleotide.PURINE, Nucleotide.R);
        Assert.assertSame(Nucleotide.PYRIMIDINE, Nucleotide.Y);
        Assert.assertSame(Nucleotide.INVALID, Nucleotide.X);
        Assert.assertSame(Nucleotide.STRONG, Nucleotide.S);
        Assert.assertSame(Nucleotide.WEAK, Nucleotide.W);
        Assert.assertSame(Nucleotide.KETO, Nucleotide.K);
        Assert.assertSame(Nucleotide.AMINO, Nucleotide.M);
    }

    @Test
    public void testLongFormNamesForTypos() {
        // We avoid a direct comparison of the constant value (rather than indirectly using its name)
        // because that would typos and here the names matter and are unlikely to change even in the long term.
        Assert.assertSame(constantNameToInstance("ADENINE"), Nucleotide.A);
        Assert.assertSame(constantNameToInstance("THYMINE"), Nucleotide.T);
        Assert.assertSame(constantNameToInstance("GUANINE"), Nucleotide.G);
        Assert.assertSame(constantNameToInstance("CYTOSINE"), Nucleotide.C);
        Assert.assertSame(constantNameToInstance("URACIL"), Nucleotide.U);
        Assert.assertSame(constantNameToInstance("ANY"), Nucleotide.N);
        Assert.assertSame(constantNameToInstance("PURINE"), Nucleotide.R);
        Assert.assertSame(constantNameToInstance("PYRIMIDINE"), Nucleotide.Y);
        Assert.assertSame(constantNameToInstance("INVALID"), Nucleotide.X);
        Assert.assertSame(constantNameToInstance("STRONG"), Nucleotide.S);
        Assert.assertSame(constantNameToInstance("WEAK"), Nucleotide.W);
        Assert.assertSame(constantNameToInstance("KETO"), Nucleotide.K);
        Assert.assertSame(constantNameToInstance("AMINO"), Nucleotide.M);
    }

    private Nucleotide constantNameToInstance(final String name) {
        try {
            return (Nucleotide) Nucleotide.class.getField(name).get(null);
        } catch (final IllegalAccessException e) {
            Assert.fail("Long name constant " + name + " is not accessible");
        } catch (final NoSuchFieldException e) {
            Assert.fail("Long name constant " + name + " does not exists");
        } catch (final ClassCastException e) {
            Assert.fail("Long name constant " + name + " so not typed as " + Nucleotide.class.getName());
        }
        throw new IllegalStateException("unreachable code");
    }

    @DataProvider(name = "testSequences")
    public Object[][] testSequences() {
        final List<Object[]> result = new ArrayList<>();
        // We add non random trivial sequences:
        result.add(new Object[]{new byte[0]});
        for (final Nucleotide nuc : Nucleotide.values()) {
            if (nuc == Nucleotide.INVALID) {
                continue;
            }
            result.add(new Object[]{new byte[]{nuc.encodeAsByte()}});
            result.add(new Object[]{Utils.repeatBytes(nuc.encodeAsByte(), MIN_RANDOM_SEQ_LENGTH)});
        }
        for (int i = 0; i < NUMBER_OF_RANDOM_SEQUENCES; i++) {
            final int length = random.nextInt(MAX_RANDOM_SEQ_LENGTH - MIN_RANDOM_SEQ_LENGTH + 1) + MIN_RANDOM_SEQ_LENGTH;
            final byte[] base = randomDNA.nextBases(length);
            result.add(new Object[]{base});
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "values")
    public Object[][] values() {
        final List<Object[]> result = new ArrayList<>(Nucleotide.values().length);
        for (final Nucleotide nuc : Nucleotide.values()) {
            result.add(new Object[]{nuc});
        }
        return result.toArray(new Object[result.size()][]);
    }
}
