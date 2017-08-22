package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class MeanQualityHistogramGeneratorUnitTest extends GATKBaseTest {

    @Test
    public void testUsingQualities() throws Exception {
        final MeanQualityByCycleSpark.HistogramGenerator hg = new MeanQualityByCycleSpark.HistogramGenerator(false);
        Assert.assertEquals(hg.useOriginalQualities, false);

        GATKRead read1 = ArtificialReadUtils.createArtificialRead("aa".getBytes(), new byte[]{50, 50}, "2M");
        hg.addRead(read1);
        assertEqualsLongArray(hg.firstReadCountsByCycle, new long[]{0, 1, 1});
        assertEqualsDoubleArray(hg.firstReadTotalsByCycle, new double[]{0, 50, 50}, 1e-05);
        assertEqualsLongArray(hg.secondReadCountsByCycle, new long[]{0, 0, 0});
        assertEqualsDoubleArray(hg.secondReadTotalsByCycle, new double[]{0, 0, 0}, 1e-05);

        GATKRead read2 = ArtificialReadUtils.createArtificialRead("aaa".getBytes(), new byte[]{11, 12, 13}, "3M");
        hg.addRead(read2);
        assertEqualsLongArray(hg.firstReadCountsByCycle, new long[]{0, 2, 2, 1});
        assertEqualsDoubleArray(hg.firstReadTotalsByCycle, new double[]{0, 61, 62, 13}, 1e-05);
        assertEqualsLongArray(hg.secondReadCountsByCycle, new long[]{0, 0, 0, 0});
        assertEqualsDoubleArray(hg.secondReadTotalsByCycle, new double[]{0, 0, 0, 0}, 1e-05);

        GATKRead read3 = ArtificialReadUtils.createArtificialRead("aa".getBytes(), new byte[]{50, 60}, "2M");
        read3.setIsReverseStrand(true);
        hg.addRead(read3);
        assertEqualsLongArray(hg.firstReadCountsByCycle, new long[]{0, 3, 3, 1});
        assertEqualsDoubleArray(hg.firstReadTotalsByCycle, new double[]{0, 121, 112, 13}, 1e-05);
        assertEqualsLongArray(hg.secondReadCountsByCycle, new long[]{0, 0, 0, 0});
        assertEqualsDoubleArray(hg.secondReadTotalsByCycle, new double[]{0, 0, 0, 0}, 1e-05);

        GATKRead read4 = ArtificialReadUtils.createArtificialRead("aaa".getBytes(), new byte[]{11, 13, 15}, "3M");
        read4.setIsReverseStrand(true);
        hg.addRead(read4);
        assertEqualsLongArray(hg.firstReadCountsByCycle, new long[]{0, 4, 4, 2});
        assertEqualsDoubleArray(hg.firstReadTotalsByCycle, new double[]{0, 121 + 15, 112 + 13, 13 + 11}, 1e-05);
        assertEqualsLongArray(hg.secondReadCountsByCycle, new long[]{0, 0, 0, 0});
        assertEqualsDoubleArray(hg.secondReadTotalsByCycle, new double[]{0, 0, 0, 0}, 1e-05);

        GATKRead read5 = ArtificialReadUtils.createArtificialRead("aaa".getBytes(), new byte[]{11, 12, 13}, "3M");
        read5.setIsSecondOfPair();
        hg.addRead(read5);
        assertEqualsLongArray(hg.firstReadCountsByCycle, new long[]{0, 4, 4, 2});
        assertEqualsDoubleArray(hg.firstReadTotalsByCycle, new double[]{0, 121 + 15, 112 + 13, 13 + 11}, 1e-05);
        assertEqualsLongArray(hg.secondReadCountsByCycle, new long[]{0, 1, 1, 1});
        assertEqualsDoubleArray(hg.secondReadTotalsByCycle, new double[]{0, 11, 12, 13}, 1e-05);

        final MeanQualityByCycleSpark.HistogramGenerator hg2 = new MeanQualityByCycleSpark.HistogramGenerator(false);
        GATKRead read1b = ArtificialReadUtils.createArtificialRead("aaaaa".getBytes(), new byte[]{51, 52, 53, 54, 55}, "5M");
        hg2.addRead(read1b);

        final MeanQualityByCycleSpark.HistogramGenerator hg2Plus1 = hg2.merge(hg);  //add short to long
        Assert.assertTrue(hg2 == hg2Plus1);
        assertEqualsLongArray(hg2Plus1.firstReadCountsByCycle, new long[]{0, 5, 5, 3, 1, 1});
        assertEqualsDoubleArray(hg2Plus1.firstReadTotalsByCycle, new double[]{0, 121 + 15 + 51, 112 + 13 + 52, 13 + 11 + 53, 54, 55}, 1e-05);
        assertEqualsLongArray(hg2Plus1.secondReadCountsByCycle, new long[]{0, 1, 1, 1, 0, 0});
        assertEqualsDoubleArray(hg2Plus1.secondReadTotalsByCycle, new double[]{0, 11, 12, 13, 0, 0}, 1e-05);

        final MeanQualityByCycleSpark.HistogramGenerator hg2_copy = new MeanQualityByCycleSpark.HistogramGenerator(false);
        hg2_copy.addRead(read1b);

        final MeanQualityByCycleSpark.HistogramGenerator hg1Plus2 = hg.merge(hg2_copy); //add long to short
        Assert.assertTrue(hg == hg1Plus2);
        assertEqualsLongArray(hg1Plus2.firstReadCountsByCycle, new long[]{0, 5, 5, 3, 1, 1});
        assertEqualsDoubleArray(hg1Plus2.firstReadTotalsByCycle, new double[]{0, 121 + 15 + 51, 112 + 13 + 52, 13 + 11 + 53, 54, 55}, 1e-05);
        assertEqualsLongArray(hg1Plus2.secondReadCountsByCycle, new long[]{0, 1, 1, 1, 0, 0});
        assertEqualsDoubleArray(hg1Plus2.secondReadTotalsByCycle, new double[]{0, 11, 12, 13, 0, 0}, 1e-05);
    }

    @Test
    public void testUsingOriginalQualities() throws Exception {
        final MeanQualityByCycleSpark.HistogramGenerator hg = new MeanQualityByCycleSpark.HistogramGenerator(true);
        Assert.assertEquals(hg.useOriginalQualities, true);

        GATKRead read1 = ArtificialReadUtils.createArtificialRead("aa".getBytes(), new byte[]{50, 50}, "2M");
        hg.addRead(read1);
        assertEqualsLongArray(hg.firstReadCountsByCycle, new long[0]);
        assertEqualsDoubleArray(hg.firstReadTotalsByCycle, new double[0], 1e-05);
        assertEqualsLongArray(hg.secondReadCountsByCycle, new long[0]);
        assertEqualsDoubleArray(hg.secondReadTotalsByCycle, new double[0], 1e-05);

        GATKRead read2 = ArtificialReadUtils.createArtificialRead("aa".getBytes(), new byte[]{50, 50}, "2M");
        read2.setAttribute(SAMTag.OQ.name(), SAMUtils.phredToFastq(new byte[]{30, 40}));
        hg.addRead(read2);

        assertEqualsLongArray(hg.firstReadCountsByCycle, new long[]{0, 1, 1});
        assertEqualsDoubleArray(hg.firstReadTotalsByCycle, new double[]{0, 30, 40}, 1e-05);
        assertEqualsLongArray(hg.secondReadCountsByCycle, new long[]{0, 0, 0});
        assertEqualsDoubleArray(hg.secondReadTotalsByCycle, new double[]{0, 0, 0}, 1e-05);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidMerge() throws Exception {
        final MeanQualityByCycleSpark.HistogramGenerator oq = new MeanQualityByCycleSpark.HistogramGenerator(true);
        final MeanQualityByCycleSpark.HistogramGenerator q = new MeanQualityByCycleSpark.HistogramGenerator(false);
        oq.merge(q);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidMergeNull() throws Exception {
        final MeanQualityByCycleSpark.HistogramGenerator oq = new MeanQualityByCycleSpark.HistogramGenerator(true);
        oq.merge(null);
    }

}