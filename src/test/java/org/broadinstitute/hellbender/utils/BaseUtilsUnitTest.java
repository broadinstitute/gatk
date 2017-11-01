package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Random;


public final class BaseUtilsUnitTest extends GATKBaseTest {
    @BeforeClass
    public void init() { }

    @Test
    public void testConvertIUPACtoN() {

        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'A', 'A'}, false, false), new byte[]{'A', 'A', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'W', 'A', 'A'}, false, false), new byte[]{'N', 'A', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'M', 'A'}, false, false), new byte[]{'A', 'N', 'A'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'A', 'A', 'K'}, false, false), new byte[]{'A', 'A', 'N'});
        checkBytesAreEqual(BaseUtils.convertIUPACtoN(new byte[]{'M', 'M', 'M'}, false, false), new byte[]{'N', 'N', 'N'});
    }

    private void checkBytesAreEqual(final byte[] b1, final byte[] b2) {
        for ( int i = 0; i < b1.length; i++ )
            Assert.assertEquals(b1[i], b2[i]);
    }

    @Test
    public void testTransitionTransversion() {
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'C', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'T', (byte)'G' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'A' ) == BaseUtils.BaseSubstitutionType.TRANSITION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'G', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );

        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'T' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'C' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'t' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'A', (byte)'c' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'t' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
        Assert.assertTrue( BaseUtils.SNPSubstitutionType( (byte)'a', (byte)'c' ) == BaseUtils.BaseSubstitutionType.TRANSVERSION );
    }

    @Test
    public void testReverseComplementString() {
        compareRCStringToExpected("ACGGT", "ACCGT");
        compareRCStringToExpected("TCGTATATCTCGCTATATATATATAGCTCTAGTATA", "TATACTAGAGCTATATATATATAGCGAGATATACGA");
        compareRCStringToExpected("AAAN", "NTTT");
    }

    private void compareRCStringToExpected(String fw, String rcExp) {
        String rcObs = new String(BaseUtils.simpleReverseComplement(fw.getBytes()));

        Assert.assertTrue(rcObs.equals(rcExp));
    }


    @Test
    public void testCaptureStdOut(){
        PrintStream stdout = System.out;
        String out = captureStdout(() -> System.out.println("Hello world"));
        GATKBaseTest.assertContains(out, "Hello world");

        try {
            captureStdout(() -> {
                throw new IllegalStateException("oh no!");
            });
        } catch (Exception e){
            //we're expecting this one... just consume it
        }
        //make sure we reset this
        Assert.assertEquals(stdout, System.out);
    }

    @Test
    public void testCaptureStderr(){
        PrintStream stderr = System.err;
        String out = captureStderr(() -> System.err.println("Hello world"));
        GATKBaseTest.assertContains(out, "Hello world");
        try {
            captureStdout(() -> {
                throw new IllegalStateException("oh no!");
            });
        } catch (Exception e){
            //we're expecting this one... just consume it
        }
        //make sure we reset this
        Assert.assertEquals(stderr, System.err);
    }

    @Test()
    public void testAssertContains(){
        assertContains("something", "thing");
        assertContains("somethingelse","some");
        assertContains("thing", "thing");

        boolean caughtException = false;
        try {
            assertContains("thing", "something"); //should fail
        } catch (AssertionError e){
            caughtException = true;
        }
        Assert.assertTrue(caughtException);
    }


    @Test(dataProvider="baseComparatorData")
    public void testBaseComparator(final Collection<byte[]> basesToSort) {
        final ArrayList<byte[]> sorted = new ArrayList<>(basesToSort);
        Collections.sort(sorted, BaseUtils.BASES_COMPARATOR);
        for (int i = 0; i < sorted.size(); i++)   {
            Assert.assertEquals(BaseUtils.BASES_COMPARATOR.compare(sorted.get(i),sorted.get(i)),0);
            final String iString = new String(sorted.get(i));
            for (int j = i; j < sorted.size(); j++) {
                final String jString = new String(sorted.get(j));
                if (iString.compareTo(jString) == 0)
                    Assert.assertEquals(BaseUtils.BASES_COMPARATOR.compare(sorted.get(i),sorted.get(j)),0);
                else
                    Assert.assertTrue(BaseUtils.BASES_COMPARATOR.compare(sorted.get(i),sorted.get(j)) * iString.compareTo(jString) > 0);
                Assert.assertTrue(BaseUtils.BASES_COMPARATOR.compare(sorted.get(i),sorted.get(j)) <= 0);
            }
        }
    }

    @DataProvider(name="baseComparatorData")
    public Object[][] baseComparatorData() {
        final int testCount = 10;
        final int testSizeAverage = 10;
        final int testSizeDeviation = 10;
        final int haplotypeSizeAverage = 100;
        final int haplotypeSizeDeviation = 100;

        final Object[][] result = new Object[testCount][];

        Utils.resetRandomGenerator();
        final Random rnd = Utils.getRandomGenerator();

        for (int i = 0; i < testCount; i++) {
            final int size = (int) Math.max(0,rnd.nextDouble() * testSizeDeviation + testSizeAverage);
            final ArrayList<byte[]> bases = new ArrayList<>(size);
            for (int j = 0; j < size; j++) {
                final int jSize = (int) Math.max(0,rnd.nextDouble() * haplotypeSizeDeviation + haplotypeSizeAverage);
                final byte[] b = new byte[jSize];
                for (int k = 0; k < jSize; k++)
                    b[k] = BaseUtils.baseIndexToSimpleBase(rnd.nextInt(4));
                bases.add(b);
            }
            result[i] = new Object[] { bases };
        }
        return result;
    }
}
