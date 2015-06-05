package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public final class CivarUnitTest extends BaseTest {


    @Test(dataProvider="validCivarExamples")
    public void testValidCivarInstanciation(final String civarString) {

        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertNotNull(civar);
    }


    @Test(dataProvider="expectedElementLengths")
    public void testValidCivarElementLength(final String civarString, final int expected) {

        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.elements().size(), expected);
    }


    @Test(dataProvider="expectedElementSizes")
    public void testValidCivarElementSizes(final String civarString, final int[] expected) {

        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.elements().size(), expected.length);
        for (int i  = 0; i < expected.length; i++) {
            Assert.assertEquals(civar.elements().get(i).size(), expected[i]);
        }
    }

    @Test(dataProvider="expectedElementOperators")
    public void testValidCivarElementOperators(final String civarString, final String expected) {

        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.elements().size(), expected.length());
        for (int i  = 0; i < expected.length(); i++) {
            Assert.assertEquals(civar.elements().get(i).operator().charValue, expected.charAt(i));
        }
    }

    @Test(dataProvider="expectedMinimumSequenceLength")
    public void testValidCivarMinimumSequenceLength(final String civarString, final int expected) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.minimumTemplateSequenceSize(), expected);
    }

    @Test(dataProvider="expectedHasVariation")
    public void testValidCivarHasVariation(final String civarString, final boolean expected) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.hasVariation(), expected);
    }


    @Test(dataProvider="invalidCivarExamples", expectedExceptions = {IllegalArgumentException.class})
    public void testInvalidInstanciation(final String civarString) {

        final Civar civar = Civar.fromCharSequence(civarString);
    }

    @Test(dataProvider="unrolledTestDataIsUnrolledExamples")
    public void testInUnrolled(final String civarString, final boolean expected) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.isUnrolled(), expected);
    }

    @Test(dataProvider="unrolledTestDataUnrolledCivarExamples")
    public void testValidCivarUnrolling(final String civarString, final String[] expected) {
        Set<String> expectedSet = new HashSet<>();
        expectedSet.addAll(Arrays.asList(expected));

        final Civar civar = Civar.fromCharSequence(civarString);
        java.util.List<Civar> unrolledList = civar.unroll();
        Assert.assertEquals(unrolledList.size(), expected.length);
        for (int i  = 0; i < expected.length; i++) {
            Assert.assertTrue(expectedSet.contains(unrolledList.get(i).toString()),
                    "Unrolled civar " + unrolledList.get(i).toString() + " not present in expected Set: " +
                            Arrays.toString(expected) + ". Unrolled set is: " + Arrays.toString(unrolledList.toArray()));
        }
    }

    @Test(dataProvider="applyToDataExamples")
    public void testValidCivarUnrolling(final String civarString, final String before, final String expectedAfter) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.applyTo(before), expectedAfter);
    }

    @Test(dataProvider="optionizeDataExamples")
    public void testValidOptionizeAll(final String civarString, final String expected) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.optionalizeAll().toString(), expected);
    }

    @DataProvider(name="validCivarExamples")
    public Iterator<Object[]> validCivarExamples() {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < VALID_CIVAR_EXAMPLES.length;
            }

            @Override
            public Object[] next() {
                return new Object[] { VALID_CIVAR_EXAMPLES[i++][0] };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="expectedHasVariation")
    public Iterator<Object[]> expectedHasVariation () {
        return validCivarExamples(5);
    }

    @DataProvider(name="expectedMinimumSequenceLength")
    public Iterator<Object[]> expectedMinimumSequenceLength () {
        return validCivarExamples(4);
    }

    @DataProvider(name="expectedElementOperators")
    public Iterator<Object[]> expectedElementOperators() {
        return validCivarExamples(3);
    }

    @DataProvider(name="expectedElementSizes")
    public Iterator<Object[]> expectedElementSizes() {
        return validCivarExamples(2);
    }

    @DataProvider(name="expectedElementLengths")
    public Iterator<Object[]> expectedElementLengths() {
        return validCivarExamples(1);
    }

    public Iterator<Object[]> validCivarExamples(final int field) {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < VALID_CIVAR_EXAMPLES.length;
            }

            @Override
            public Object[] next() {
                return new Object[] { VALID_CIVAR_EXAMPLES[i][0], VALID_CIVAR_EXAMPLES[i++][field] };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="unrolledTestDataIsUnrolledExamples")
    public Iterator<Object[]> unrolledTestDataIsUnrolledExamples() {
        return unrolledTestDataExamples(1);
    }

    @DataProvider(name="unrolledTestDataUnrolledCivarExamples")
    public Iterator<Object[]> unrolledTestDataUnrolledCivarExamples() {
        return unrolledTestDataExamples(2);
    }

    public Iterator<Object[]> unrolledTestDataExamples(final int field) {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < UNROLLED_TEST_DATA.length;
            }

            @Override
            public Object[] next() {
                return new Object[] { UNROLLED_TEST_DATA[i][0], UNROLLED_TEST_DATA[i++][field] };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="optionizeDataExamples")
    public Iterator<Object[]> optionizeDataExamples() {
        return optionizeDataExamples(1);
    }

    public Iterator<Object[]> optionizeDataExamples(final int field) {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < OPTIONIZED_TEST_DATA.length;
            }

            @Override
            public Object[] next() {
                return new Object[] { OPTIONIZED_TEST_DATA[i][0], OPTIONIZED_TEST_DATA[i++][field] };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="applyToDataExamples")
    public Iterator<Object[]> applyToDataExamples() {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < APPLY_TO_TEST_DATA.length;
            }

            @Override
            public Object[] next() {
                return APPLY_TO_TEST_DATA[i++];
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="invalidCivarExamples")
    public Object[][] invalidCivarExamples() {
            return INVALID_CIVAR_EXAMPLES;
    }

    // columns : Civar string, number of elements.
    private static final Object[][] INVALID_CIVAR_EXAMPLES = new Object[][] {
            {"(100="},
            {"*=)"},
            {"10(=2T30="},
            {"2*=2T/3*="},
            {"3I(acc)"},
            {"a"},
            {")"},
            {"100&=1"},
            {"?100="},

    };


    private static final Object[][] VALID_CIVAR_EXAMPLES = new Object[][] {
            {"100=", 1, ints(100), "=", 100, false },
            {"*=", 1 , ints(1), "=", 0, false },
            {"10=2T30=", 3, ints(10,2,30), "=T=",42 , true},
            {"*=2T3*=", 3, ints(1,2,3), "=T=",2 , true},
            {"3Iacc",1 , ints(3), "I", 0, true},
            {"Ia",1, ints(1), "I", 0, true},
            {"10D",1, ints(10), "D", 10, true},
            {"*", 1, ints(1), "=", 0, false},
            {"*D", 1, ints(1), "D", 0, true},
            {"10(1D)10=",3, ints(10,1,10), "=(=", 21, true},
            {"1*",1, ints(1), "=", 0, false},
            {"1*2*",2, ints(1,2), "==", 0, false},
            {"*11",2, ints(1,11), "==", 11, false},
            {"100=1T100=", 3, ints(100,1,100), "=T=", 201, true},
            {"100=3Iacg101=", 3, ints(100,3,101), "=I=", 201, true},
            {"100=30Igctcggatgccttgcggggctccagagtcc101=", 3 , ints(100,30,101), "=I=", 201, true},
            {"99=3D99=", 3, ints(99,3,99), "=D=", 201, true},
            {"84=30D84=", 3, ints(84,30,84), "=D=", 198, true},
            {"91=1T9=3Iacg100=", 5, ints(91,1,9,3,100), "=T=I=", 201, true},
            {"71=1T29=3Iacg100=",5, ints(71,1,29,3,100), "=T=I=",201, true},
            {"75=1T8=1T8=1T8=1T8=1T75=", 11, ints(75,1,8,1,8,1,8,1,8,1,75), "=T=T=T=T=T=",187, true},
            {"75=1T?8=", 3, ints(75,1,8), "=T=", 84, true}
    };

    private static final Object[][] UNROLLED_TEST_DATA = new Object[][] {
            { "10=1D10=", true, strs( "10=1D10=") },
            { "10=(1D)10=", true, strs( "10=(1D)10=") },
            { "10=1D?10=", false, strs("10=1=10=", "10=1D10=") },
            { "10=1D?10=3Iacg?10=", false , strs("10=1=10=0=10=","10=1=10=3Iacg10=", "10=1D10=0=10=", "10=1D10=3Iacg10=") },
            {  "10=1D?10=" , false, strs("10=1D10=","10=1=10=") },
            {  "100=1T?100=" , false, strs("100=1T100=","100=1=100=") },
            {  "100=3Iacg?101=" , false, strs("100=3Iacg101=","100=0=101=") },
            {  "100=30Igctcggatgccttgcggggctccagagtcc?101=", false ,strs("100=30Igctcggatgccttgcggggctccagagtcc101=", "100=0=101=") },
            {  "99=3D?99=", false , strs("99=3D99=","99=3=99=") },
            {  "84=30D?84=", false, strs("84=30D84=", "84=30=84=")},
            {  "91=1T?9=3Iacg?100=", false, strs("91=1T9=3Iacg100=", "91=1=9=3Iacg100=", "91=1=9=0=100=", "91=1T9=0=100=") },
            {  "71=1T?29=3Iacg?100=", false , strs("71=1T29=3Iacg100=","71=1=29=3Iacg100=","71=1=29=0=100=", "71=1T29=0=100=") },
           // {  "75=1T?8=1T?8=1T?8=1T?8=1T?75=", false, },
            {  "75=1T?8=", false, strs("75=1T8=","75=1=8=") }
    };

    private static final Object[][] OPTIONIZED_TEST_DATA = new Object[][] {
            { "10=1D10=", "10=1D?10=" },
            {"100=1T100=","100=1T?100=" },
            {"100=3Iacg101=", "100=3Iacg?101=" },
            {"100=30Igctcggatgccttgcggggctccagagtcc101=","100=30Igctcggatgccttgcggggctccagagtcc?101="},
            {"99=3D99=", "99=3D?99="},
            {"84=30D84=", "84=30D?84="},
            {"91=1T9=3Iacg100=", "91=1T?9=3Iacg?100="},
            {"71=1T29=3Iacg100=","71=1T?29=3Iacg?100="},
            {"75=1T8=1T8=1T8=1T8=1T75=", "75=1T?8=1T?8=1T?8=1T?8=1T?75="},
            {"75=1T?8=", "75=1T?8="}
    };

    private static final Object[][] APPLY_TO_TEST_DATA = new Object[][] {
            {"3=1D3=", "ACTAACT", "ACTACT" },
            {"*=1C*=","ACTTACT", "ACTAACT" },
            {"4=3Iacg3=","ACTGACT","ACTGACGACT" },
            {"*=30Igctcggatgccttgcggggctccagagtcc*=","AA","AGCTCGGATGCCTTGCGGGGCTCCAGAGTCCA"},
            {"*=3D*=", "ACTTTAC","ACAC"},
            {"1=30D1=", "AGCTCGGATGCCTTGCGGGGCTCCAGAGTCCA","AA"},
    };


    private static int[] ints(final int ... iii) {
        return iii;
    }

    private static String[] strs(final String... sss) {
        return sss;
    }

}
