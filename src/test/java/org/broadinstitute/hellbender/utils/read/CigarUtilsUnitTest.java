package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import org.broadinstitute.hellbender.testutils.ReadClipperTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public final class CigarUtilsUnitTest {

    @DataProvider(name = "testData_invertCigar")
    public Iterator<Object[]> testData_invertCigar(final Method testMethod) {
        final String[][] TEST_CIGARS = {
                {"10M", "10M"},
                {"1M2M", "2M1M"},
                {"10M10D", "10D10M"},
                {"2I2D", "2D2I"},
                {"2D2I", "2I2D"},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_invertCigar")
    public void testInvertCigar(final String cigarStrIn, final String expectedCigarStrOut){
        final Cigar cigarIn = TextCigarCodec.decode(cigarStrIn);
        final Cigar cigarOut = CigarUtils.invertCigar(cigarIn);
        final String actualCigarStrOut = TextCigarCodec.encode(cigarOut);
        Assert.assertEquals(actualCigarStrOut, expectedCigarStrOut);
    }

    @DataProvider(name = "testData_unclipCigar")
    public Iterator<Object[]> testData_unclipCigar(final Method testMethod) {
        final String[][] TEST_CIGARS = {
                {"10M", "10M"},
                {"1M1D1M1I1M", "1M1D1M1I1M"},
                {"1H10M", "10M"},
                {"1H1S10M", "10M"},
                {"10M1S1H", "10M"},
                {"1H1S10M1S1H", "10M"},
                {"1H1S10M1D10M1S1H", "10M1D10M"},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_unclipCigar")
    public void testUnclipCigar(final String cigarStrIn, final String expectedCigarStrOut){
        final Cigar cigarIn = TextCigarCodec.decode(cigarStrIn);
        final Cigar cigarOut = CigarUtils.removeClipsAndPadding(cigarIn);
        final String actualCigarStrOut = TextCigarCodec.encode(cigarOut);
        Assert.assertEquals(actualCigarStrOut, expectedCigarStrOut);
    }

    @DataProvider(name = "testData_countRefBasesBasedOnCigar")
    public Iterator<Object[]> testData_countRefBasesBasedOnCigar(final Method testMethod) {
        final Object[][] TEST_CIGARS = {
                {"10M", 0, 1, 10},
                {"10M1D", 0, 1, 10},
                {"10M1D1S", 0, 1, 10},
                {"10M1D1N1S", 0, 1, 10},
                {"10M1D1N1S1H", 0, 1, 10},
                {"10M1D", 0, 2, 11},
                {"10M1D1S", 0, 2, 11},
                {"10M1D1N1S", 0, 2, 11},
                {"10M1D1N1S1H", 0, 2, 11},
                {"10M1=1X1S", 0, 2, 11},
                {"10M1=1X1S1H", 0, 2, 11},
                {"10M1D2N4S", 1, 3, 3},
                {"10M1D2N4S8H", 1, 4, 7},
                {"10M1I2N4S", 1, 3, 2},
                {"10M1I2N4S8H", 1, 4, 6},

                //https://github.com/broadinstitute/hellbender/issues/450
                {"1M1S1H1H", 1, 2, 1},
                {"1M1I1H1H", 1, 2, 0},//I is skipped
                {"1M1P1H1H", 1, 2, 0},//P is skipped
                {"1M1H1H1H", 1, 2, 1},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_countRefBasesBasedOnCigar")
    public void testCountRefBasesBasedOnCigar(final String cigarStrIn, final int start, final int end, final int expected){
        final int actual = CigarUtils.countRefBasesAndClips(TextCigarCodec.decode(cigarStrIn).getCigarElements(), start, end);
        Assert.assertEquals(actual, expected, cigarStrIn);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCountRefBasesBasedOnCigarNull(){
        CigarUtils.countRefBasesAndClips(null, 1, 2);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCountRefBasesBasedOnCigarStart1(){
        final String cigarStrIn = "1M1=1X";
        CigarUtils.countRefBasesAndClips(TextCigarCodec.decode(cigarStrIn).getCigarElements(), -1, 1);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCountRefBasesBasedOnCigarStart2(){
        final String cigarStrIn = "1M1=1X";
        CigarUtils.countRefBasesAndClips(TextCigarCodec.decode(cigarStrIn).getCigarElements(), 2, 1);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCountRefBasesBasedOnCigarEnd2(){
        final String cigarStrIn = "1M1=1X";
        CigarUtils.countRefBasesAndClips(TextCigarCodec.decode(cigarStrIn).getCigarElements(), 1, 6);
    }


    @DataProvider(name = "testData_testComputeCigar")
    public Iterator<Object[]> testData_testComputeCigar(final Method testMethod) {
        final String[][] TEST_CIGARS = {
                {"ATGGAGGGGC", "ATGGTGGGGC", "10M"},
                {"ATGGAGGGGC", "ATGGAAAATGGGGC", "5M4I5M"},
                {"ATGGAGGGGC", "ATGGAAAAAAAAATGGGGC", "5M9I5M"},
                {"ATGGAAAAAGGGGC", "ATGGTGGGGC", "4M4D6M"},
                {"ATGGAAAAAGGGGC", "ATGGAAAATGGGGC", "14M"},
                {"ATGGAAAAAGGGGC", "ATGGAAAAAAAAATGGGGC", "9M5I5M"},
                {"ATGGAAAAAAAAAAGGGGC", "ATGGTGGGGC", "4M9D6M"},
                {"ATGGAAAAAAAAAAGGGGC", "ATGGAAAATGGGGC", "4M5D10M"},
                {"ATGGAAAAAAAAAAGGGGC", "ATGGAAAAAAAAATGGGGC", "19M"},
                {"NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN", "NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN", "51M"},
                {"NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN", "NNNACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN", "3M6I48M"},
                {"NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN", "NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN", "51M"},
                {"NNNTGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN", "NNNACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGANNN", "3M6I48M"},
                {"TCCCCCGGGT", "TCCCCCGGGT", "10M"},
                {"TCCCCCGGGT", "TAAACCCCCT", "1M3I5M3D1M"},
                {"G", "C", "1M"},
                {"G", "", "1D"},
                {"", "C", "1I"},
                {"AAA", "CGT", "3M"},
                {"TAT", "CAC", "3M"},
                {"GCTG", "GTCG", "4M"},
                {"AAAAA", "", "5D"},
                {"", "AAAAA", "5I"},
                {"AAAAACC", "CCGGGGGG", "5D2M6I"},
                {"GX", "CX", "2M"},
                {"GX", "X", "1D1M"},
                {"X", "CX", "1I1M"},
                {"AAAX", "CGTX", "4M"},
                {"TATX", "CACX", "4M"},
                {"GCTGX", "GTCGX", "5M"},
                {"AAAAAX", "X", "5D1M"},
                {"X", "AAAAAX", "5I1M"},
                {"AAAAACCX", "CCGGGGGGX", "5D2M6I1M"},
                {"GXXXXXXXXXXXXX", "CXXXXXXXXXXXXX", "14M"},
                {"GXXXXXXXXXXXXX", "XXXXXXXXXXXXX", "1D13M"},
                {"XXXXXXXXXXXXX", "CXXXXXXXXXXXXX", "1I13M"},
                {"AAAXXXXXXXXXXXXX", "CGTXXXXXXXXXXXXX", "16M"},
                {"TATXXXXXXXXXXXXX", "CACXXXXXXXXXXXXX", "16M"},
                {"GCTGXXXXXXXXXXXXX", "GTCGXXXXXXXXXXXXX", "17M"},
                {"AAAAAXXXXXXXXXXXXX", "XXXXXXXXXXXXX", "5D13M"},
                {"XXXXXXXXXXXXX", "AAAAAXXXXXXXXXXXXX", "5I13M"},
                {"AAAAACCXXXXXXXXXXXXX", "CCGGGGGGXXXXXXXXXXXXX", "5D2M6I13M"},
                {"XG", "XC", "2M"},
                {"XG", "X", "1M1D"},
                {"X", "XC", "1M1I"},
                {"XAAA", "XCGT", "4M"},
                {"XTAT", "XCAC", "4M"},
                {"XGCTG", "XGTCG", "5M"},
                {"XAAAAA", "X", "1M5D"},
                {"X", "XAAAAA", "1M5I"},
                {"XAAAAACC", "XCCGGGGGG", "1M5D2M6I"},
                {"XGX", "XCX", "3M"},
                {"XGX", "XX", "1M1D1M"},
                {"XX", "XCX", "1M1I1M"},
                {"XAAAX", "XCGTX", "5M"},
                {"XTATX", "XCACX", "5M"},
                {"XGCTGX", "XGTCGX", "6M"},
                {"XAAAAAX", "XX", "1M5D1M"},
                {"XX", "XAAAAAX", "1M5I1M"},
                {"XAAAAACCX", "XCCGGGGGGX", "1M5D2M6I1M"},
                {"XGXXXXXXXXXXXXX", "XCXXXXXXXXXXXXX", "15M"},
                {"XGXXXXXXXXXXXXX", "XXXXXXXXXXXXXX", "1M1D13M"},
                {"XXXXXXXXXXXXXX", "XCXXXXXXXXXXXXX", "1M1I13M"},
                {"XAAAXXXXXXXXXXXXX", "XCGTXXXXXXXXXXXXX", "17M"},
                {"XTATXXXXXXXXXXXXX", "XCACXXXXXXXXXXXXX", "17M"},
                {"XGCTGXXXXXXXXXXXXX", "XGTCGXXXXXXXXXXXXX", "18M"},
                {"XAAAAAXXXXXXXXXXXXX", "XXXXXXXXXXXXXX", "1M5D13M"},
                {"XXXXXXXXXXXXXX", "XAAAAAXXXXXXXXXXXXX", "1M5I13M"},
                {"XAAAAACCXXXXXXXXXXXXX", "XCCGGGGGGXXXXXXXXXXXXX", "1M5D2M6I13M"},
                {"XXXXXXXXXXXXXG", "XXXXXXXXXXXXXC", "14M"},
                {"XXXXXXXXXXXXXG", "XXXXXXXXXXXXX", "13M1D"},
                {"XXXXXXXXXXXXX", "XXXXXXXXXXXXXC", "13M1I"},
                {"XXXXXXXXXXXXXAAA", "XXXXXXXXXXXXXCGT", "16M"},
                {"XXXXXXXXXXXXXTAT", "XXXXXXXXXXXXXCAC", "16M"},
                {"XXXXXXXXXXXXXGCTG", "XXXXXXXXXXXXXGTCG", "17M"},
                {"XXXXXXXXXXXXXAAAAA", "XXXXXXXXXXXXX", "13M5D"},
                {"XXXXXXXXXXXXX", "XXXXXXXXXXXXXAAAAA", "13M5I"},
                {"XXXXXXXXXXXXXAAAAACC", "XXXXXXXXXXXXXCCGGGGGG", "13M5D2M6I"},
                {"XXXXXXXXXXXXXGX", "XXXXXXXXXXXXXCX", "15M"},
                {"XXXXXXXXXXXXXGX", "XXXXXXXXXXXXXX", "13M1D1M"},
                {"XXXXXXXXXXXXXX", "XXXXXXXXXXXXXCX", "13M1I1M"},
                {"XXXXXXXXXXXXXAAAX", "XXXXXXXXXXXXXCGTX", "17M"},
                {"XXXXXXXXXXXXXTATX", "XXXXXXXXXXXXXCACX", "17M"},
                {"XXXXXXXXXXXXXGCTGX", "XXXXXXXXXXXXXGTCGX", "18M"},
                {"XXXXXXXXXXXXXAAAAAX", "XXXXXXXXXXXXXX", "13M5D1M"},
                {"XXXXXXXXXXXXXX", "XXXXXXXXXXXXXAAAAAX", "13M5I1M"},
                {"XXXXXXXXXXXXXAAAAACCX", "XXXXXXXXXXXXXCCGGGGGGX", "13M5D2M6I1M"},
                {"XXXXXXXXXXXXXGXXXXXXXXXXXXX", "XXXXXXXXXXXXXCXXXXXXXXXXXXX", "27M"},
                {"XXXXXXXXXXXXXGXXXXXXXXXXXXX", "XXXXXXXXXXXXXXXXXXXXXXXXXX", "13M1D13M"},
                {"XXXXXXXXXXXXXXXXXXXXXXXXXX", "XXXXXXXXXXXXXCXXXXXXXXXXXXX", "13M1I13M"},
                {"XXXXXXXXXXXXXAAAXXXXXXXXXXXXX", "XXXXXXXXXXXXXCGTXXXXXXXXXXXXX", "29M"},
                {"XXXXXXXXXXXXXTATXXXXXXXXXXXXX", "XXXXXXXXXXXXXCACXXXXXXXXXXXXX", "29M"},
                {"XXXXXXXXXXXXXGCTGXXXXXXXXXXXXX", "XXXXXXXXXXXXXGTCGXXXXXXXXXXXXX", "30M"},
                {"XXXXXXXXXXXXXAAAAAXXXXXXXXXXXXX", "XXXXXXXXXXXXXXXXXXXXXXXXXX", "13M5D13M"},
                {"XXXXXXXXXXXXXXXXXXXXXXXXXX", "XXXXXXXXXXXXXAAAAAXXXXXXXXXXXXX", "13M5I13M"},
                {"XXXXXXXXXXXXXAAAAACCXXXXXXXXXXXXX", "XXXXXXXXXXXXXCCGGGGGGXXXXXXXXXXXXX", "13M5D2M6I13M"},
                {"ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG", "ATGGTTTCCTCCCCAAAAAAAAAAAATGGCCGCCCC", "32M4I"},
                {"ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG", "ATGGTTTCCTCCCCAAAAAAAAAAAATGGCCG", "32M"},
                {"ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG", "ATGGAAAATTTCCTCCCCAAAAAAAAAAAAGGGGTGGCCGCCCC", "5M4I22M4I5M4I"},
                {"ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG", "ATGGAAAATTTCCTCCCCAAAAAAAAAAAAGGGGTGGCCG", "5M4I22M4I5M"},
                {"ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG", "ATGGAAAAAAAAATTTCCTCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCGCCCC", "5M9I22M9I5M4I"},
                {"ATGGATTCCTCCCCAAAAAAAAAAAAGGGCCG", "ATGGAAAAAAAAATTTCCTCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCG", "5M9I22M9I5M"},
                {"ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG", "ATGGTTTCCTCCCCCCCCAAAAAAAAAAAATGGCCGCCCC", "4M4D26M4D6M4I"},
                {"ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG", "ATGGTTTCCTCCCCCCCCAAAAAAAAAAAATGGCCG", "4M4D26M4D6M"},
                {"ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG", "ATGGAAAATTTCCTCCCCCCCCAAAAAAAAAAAAGGGGTGGCCGCCCC", "44M4I"},
                {"ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG", "ATGGAAAATTTCCTCCCCCCCCAAAAAAAAAAAAGGGGTGGCCG", "44M"},
                {"ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG", "ATGGAAAAAAAAATTTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCGCCCC", "9M5I30M5I5M4I"},
                {"ATGGAAAAATTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGCCG", "ATGGAAAAAAAAATTTCCTCCCCCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCG", "9M5I30M5I5M"},
                {"ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG", "ATGGTTTCCTCCCCCCCCCCCCCAAAAAAAAAAAATGGCCGCCCC", "4M9D31M9D6M4I"},
                {"ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG", "ATGGTTTCCTCCCCCCCCCCCCCAAAAAAAAAAAATGGCCG", "4M9D31M9D6M"},
                {"ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG", "ATGGAAAATTTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGTGGCCGCCCC", "4M5D35M5D10M4I"},
                {"ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG", "ATGGAAAATTTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGTGGCCG", "4M5D35M5D10M"},
                {"ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG", "ATGGAAAAAAAAATTTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCGCCCC", "59M4I"},
                {"ATGGAAAAAAAAAATTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGGGGCCG", "ATGGAAAAAAAAATTTCCTCCCCCCCCCCCCCAAAAAAAAAAAAGGGGGGGGGTGGCCG", "59M"},
        };
        final List<Object[]> result = new LinkedList<>();
        Collections.addAll(result, TEST_CIGARS);
        return result.iterator();
    }

    @Test(dataProvider = "testData_testComputeCigar")
    public void testComputeCigar(String s1, String s2, String expectedCigar) throws Exception {
        final Cigar actualCigar = CigarUtils.calculateCigar(s1.getBytes(), s2.getBytes(), SmithWatermanJavaAligner.getInstance(), SWOverhangStrategy.SOFTCLIP);
        final Cigar decode = TextCigarCodec.decode(expectedCigar);
        Assert.assertEquals(actualCigar, decode);
    }

    @Test(dataProvider = "randomValidCigars")
    public void testLeftClip(final Cigar cigar) {
        final int actual = CigarUtils.countLeftClippedBases(cigar);
        int expected = 0;
        for (final CigarElement element : cigar.getCigarElements()) {
            if (!element.getOperator().isClipping()) {
                break;
            }
            expected += element.getLength();
        }
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "allClipFunkyCigars", expectedExceptions = IllegalArgumentException.class)
    public void testLeftClipFunkly(final Cigar cigar) {
        CigarUtils.countLeftClippedBases(cigar);
    }

    @Test
    public void testLeftClipByHand() {
        final List<Pair<Cigar, Integer>> tests =  Stream.of(
                "13H3M35D13M2I45S30H, 13",
                "1H3S67M13S, 4",
                "13M30S, 0",
                "113S4M12S, 113")
                .map(s -> s.split(",\\s*"))
                .map(ss -> new ImmutablePair<>(TextCigarCodec.decode(ss[0]), Integer.parseInt(ss[1]) ))
                .collect(Collectors.toList());
        for (final Pair<Cigar, Integer> test : tests) {
            Assert.assertEquals(CigarUtils.countLeftClippedBases(test.getLeft()), test.getRight().intValue(), "" + test.getLeft());
        }
    }

    @Test(dataProvider = "randomValidCigars")
    public void testLeftHardClip(final Cigar cigar) {
        final int actual = CigarUtils.countLeftHardClippedBases(cigar);
        int expected = 0;
        for (final CigarElement element : cigar.getCigarElements()) {
            if (element.getOperator() != CigarOperator.H) {
                break;
            }
            expected += element.getLength();
        }
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testLeftHardClipByHand() {
        final List<Pair<Cigar, Integer>> tests =  Stream.of(
                "13H3M35D13M2I45S30H, 13",
                "1H3S67M13S, 1",
                "13M30S, 0",
                "113S4M12S, 0")
                .map(s -> s.split(",\\s*"))
                .map(ss -> new ImmutablePair<>(TextCigarCodec.decode(ss[0]), Integer.parseInt(ss[1]) ))
                .collect(Collectors.toList());
        for (final Pair<Cigar, Integer> test : tests) {
            Assert.assertEquals(CigarUtils.countLeftHardClippedBases(test.getLeft()), test.getRight().intValue(), "" + test.getLeft());
        }
    }

    @Test(dataProvider = "randomValidCigars")
    public void testRightClip(final Cigar cigar) {
        final int actual = CigarUtils.countRightClippedBases(cigar);
        int expected = 0;
        boolean nonClippingFound = false;
        for (final CigarElement element : CigarUtils.invertCigar(cigar).getCigarElements()) {
            if (!element.getOperator().isClipping()) {
                nonClippingFound = true;
                break;
            }
            expected += element.getLength();
        }
        Assert.assertEquals(actual, nonClippingFound ? expected : 0);
    }

    @Test(dataProvider = "allClipFunkyCigars", expectedExceptions = IllegalArgumentException.class)
    public void testRightClipFunkly(final Cigar cigar) {
        CigarUtils.countRightClippedBases(cigar);
    }

    @Test(dataProvider = "randomValidCigars")
    public void testRightHardClip(final Cigar cigar) {
        final int actual = CigarUtils.countRightHardClippedBases(cigar);
        int expected = 0;
        for (final CigarElement element : CigarUtils.invertCigar(cigar).getCigarElements()) {
            if (element.getOperator() != CigarOperator.H) {
                break;
            }
            expected += element.getLength();
        }
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testRightHardClipByHand() {
        final List<Pair<Cigar, Integer>> tests =  Stream.of(
                "13H3M35D13M2I45S30H, 30",
                "1H3S67M13S, 0",
                "13M30S10S1H, 1",
                "113S4M, 0")
                .map(s -> s.split(",\\s*"))
                .map(ss -> new ImmutablePair<>(TextCigarCodec.decode(ss[0]), Integer.parseInt(ss[1]) ))
                .collect(Collectors.toList());
        for (final Pair<Cigar, Integer> test : tests) {
            Assert.assertEquals(CigarUtils.countRightHardClippedBases(test.getLeft()), test.getRight().intValue(), "" + test.getLeft());
        }
    }

    @Test(dataProvider = "randomValidCigars")
    public void testReadLength(final Cigar cigar) {
        final int actual = CigarUtils.countUnclippedReadBases(cigar);
        final int expected = cigar.getCigarElements().stream()
                .filter(ce -> ce.getOperator().consumesReadBases() || ce.getOperator().isClipping())
                .mapToInt(CigarElement::getLength).sum();
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name = "allClipFunkyCigars")
    public Object[][] allClipFunkyCigars() {
        final Object[][] randomValidCigars = randomValidCigars();
        final List<Object[]> result = new ArrayList<>(randomValidCigars.length + 10);
        result.add(new Object[] {TextCigarCodec.decode("10S20S30S")}); // all soft-clips.
        result.add(new Object[] {TextCigarCodec.decode("10H10S")}); // H and S.
        result.add(new Object[] {TextCigarCodec.decode("10H10S30H")}); // H, S and H.
        result.add(new Object[] {TextCigarCodec.decode("10S30H")}); // S and H.
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name = "randomValidCigars")
    public static Object[][] randomValidCigars() {
        return CigarTestUtils.randomValidCigars(new Random(13), 1_000,
                10, 100, new Cigar()).stream()
                .map(x -> new Object[] { x }).toArray(Object[][]::new);
    }

    @Test
    public void testCountAlignedBases() {
        Assert.assertEquals(CigarUtils.countAlignedBases(TextCigarCodec.decode("548M1215H")), 548);
        Assert.assertEquals(CigarUtils.countAlignedBases(TextCigarCodec.decode("545H252M966H")), 252);
        Assert.assertEquals(CigarUtils.countAlignedBases(TextCigarCodec.decode("797S966M")), 966);
        Assert.assertEquals(CigarUtils.countAlignedBases(TextCigarCodec.decode("79H113M14D15M16D16M3I17M2I117M797H")), 278);
    }

    @DataProvider(name = "revert_soft_clips")
    public Object[][] revertSoftClips() {
        return new Object[][] {
                {"10M","10M"},
                {"10H10M", "10H10M"},
                {"10S10M", "20M"},
                {"10S10M10S", "30M"},
                {"10S10M10S10H", "30M10H"}
        };
    }

    @Test(dataProvider = "revert_soft_clips")
    public void testRevertSoftClips(final String original, final String expected) {
        Assert.assertEquals(CigarUtils.revertSoftClips(TextCigarCodec.decode(original)).toString(), expected);
    }

    @DataProvider(name = "clip_cigar")
    public Object[][] ClipCigar() {
        return new Object[][] {
                //simple cases
                {"10M", 0, 5, "5S5M", "5H5M"},
                {"10M", 5, 10, "5M5S", "5M5H"},
                {"10H10M", 0, 5, "10H5S5M", "15H5M"},
                {"10H10M", 5, 10, "10H5M5S", "10H5M5H"},
                {"10M10H", 0, 5, "5S5M10H", "5H5M10H"},

                // clipping into insertion
                {"10M10I10M", 0, 5, "5S5M10I10M", "5H5M10I10M"},
                {"10M10I10M", 0, 15, "15S5I10M", "15H5I10M"},
                {"10M10I10M", 15, 30, "10M5I15S", "10M5I15H"},

                // clipping into a soft clip
                {"10S10M10S", 0, 5, "10S10M10S", "5H5S10M10S"},
                {"10S10M10S", 25, 30, "10S10M10S", "10S10M5S5H"},
                {"10S10M10S", 0, 15, "15S5M10S", "15H5M10S"},

                // clipping over a deletion
                {"10M10D10M", 0, 10, "10S10M", "10H10M"},
                {"10M10D10M", 0, 15, "15S5M", "15H5M"},
                {"10M10D10M", 5, 20, "5M15S", "5M15H"},

                // removing leading deletions
                {"10D10M", 0, 5, "5S5M", "5H5M"}
        };
    }

    @Test(dataProvider = "clip_cigar")
    public void testClipCigar(final String original, final int start, final int stop, final String expectedSoftClip, final String expectedHardClip) {
        Assert.assertEquals(CigarUtils.clipCigar(TextCigarCodec.decode(original), start, stop, CigarOperator.SOFT_CLIP).toString(), expectedSoftClip);
        Assert.assertEquals(CigarUtils.clipCigar(TextCigarCodec.decode(original), start, stop, CigarOperator.HARD_CLIP).toString(), expectedHardClip);
    }

    @DataProvider(name = "alignment_start_shift")
    public Object[][] alignmentStartShift() {
        return new Object[][] {
                {"70M", 10, 10},
                {"70M", 0, 0},

                {"30M10D30M", 29, 29},
                {"30M10D30M", 30, 40},
                {"30M10D30M", 31, 41},

                {"30M10D30M", 29, 29},
                {"30M10I30M", 30, 30},
                {"30M10I30M", 31, 30},
                {"30M10I30M", 40, 30},
                {"30M10I30M", 41, 31},

                {"10H10M", 5, 5},
                {"10S10M", 5, 0},
                {"10S10M", 5, 0},
        };
    }

    @Test(dataProvider = "alignment_start_shift")
    public void testAlignmentStartShift(final String cigarString, final int numClips, final int expectedResult) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final int actualResult = CigarUtils.alignmentStartShift(cigar, numClips);
        Assert.assertEquals(actualResult, expectedResult);
    }
}
