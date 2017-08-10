package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;
import scala.Tuple3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class VariantDetectorFromLongReadAlignmentsForSimpleStrandSwitchUnitTest extends BaseTest {

    @DataProvider(name = "forComputeNewRefSpanAndCigar")
    private Object[][] createTestDataForComputeNewRefSpanAndCigar() throws IOException {

        final List<Object[]> data = new ArrayList<>(20);

        AlignmentInterval alignment = new AlignmentInterval(new SimpleInterval("chr1", 175417007, 175417074),
                14, 81, TextCigarCodec.decode("13H68M394H"),
                true, 60, 0, 68, false, false);
        SimpleInterval refSpan = new SimpleInterval("chr1", 175417007, 175417043);
        data.add(new Object[]{alignment, 31, true, refSpan, TextCigarCodec.decode("13H37M31S394H")});

        alignment = new AlignmentInterval(new SimpleInterval("chr2", 122612655, 122612751),
                9, 105, TextCigarCodec.decode("8H97M138H"),
                false, 60, 0, 97, false, false);
        refSpan = new SimpleInterval("chr2", 122612659, 122612751);
        data.add(new Object[]{alignment, 4, true, refSpan, TextCigarCodec.decode("8H93M4S138H")});

        alignment = new AlignmentInterval(new SimpleInterval("chr6", 66782514, 66782679),
                32, 197, TextCigarCodec.decode("31S166M"),
                false, 60, 3, 151, false, false);
        refSpan = new SimpleInterval("chr6", 66782514, 66782675);
        data.add(new Object[]{alignment, 4, false, refSpan, TextCigarCodec.decode("35S162M")});

        alignment = new AlignmentInterval(new SimpleInterval("chr2", 91421528, 91421734),
                271, 477, TextCigarCodec.decode("270H207M"),
                true, 40, 12, 147, false, false);
        refSpan = new SimpleInterval("chr2", 91421560, 91421734);
        data.add(new Object[]{alignment, 32, false, refSpan, TextCigarCodec.decode("270H32S175M")});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "forComputeNewRefSpanAndCigar", groups = "sv")
    public void testComputeNewRefSpanAndCigar(final AlignmentInterval interval, final int clipLength, final boolean clipFrom3PrimeEnd,
                                              final SimpleInterval expectedRefSpan, final Cigar expectedCigar) {

        final Tuple2<SimpleInterval, Cigar> x = SimpleStrandSwitchVariantDetector.computeNewRefSpanAndCigar(interval, clipLength, clipFrom3PrimeEnd);
        Assert.assertEquals(x._1, expectedRefSpan);
        Assert.assertEquals(x._2, expectedCigar);
    }

    @DataProvider(name = "forCigarExtraction")
    private Object[][] createTestDataForCigarExtraction() throws IOException {

        final List<Object[]> data = new ArrayList<>(20);

        SimpleInterval refSpan = new SimpleInterval("chr1", 82666357, 82666765);
        AlignmentInterval alignment = new AlignmentInterval(refSpan, 69, 472,
                TextCigarCodec.decode("68S122M5D282M"), true, 60, 11,
                353, false, false);
        List<CigarElement> left = Arrays.asList(new CigarElement(68, CigarOperator.S));
        List<CigarElement> middle = Arrays.asList(new CigarElement(122, CigarOperator.M), new CigarElement(5, CigarOperator.D), new CigarElement(282, CigarOperator.M));
        List<CigarElement> right = Collections.emptyList();
        data.add(new Object[]{alignment, left, middle, right});

        refSpan = new SimpleInterval("chr3", 61792401, 61792448);
        alignment = new AlignmentInterval(refSpan, 43, 90,
                TextCigarCodec.decode("42H48M382H"), true, 46, 1,
                43, false, false);
        left = Arrays.asList(new CigarElement(42, CigarOperator.H));
        middle = Arrays.asList(new CigarElement(48, CigarOperator.M));
        right = Arrays.asList(new CigarElement(382, CigarOperator.H));
        data.add(new Object[]{alignment, left, middle, right});

        refSpan = new SimpleInterval("chrY", 26303624, 26303671);
        alignment = new AlignmentInterval(refSpan, 1, 48,
                TextCigarCodec.decode("48M424H"), true, 0, 2,
                38, false, false);
        left = Collections.emptyList();
        middle = Arrays.asList(new CigarElement(48, CigarOperator.M));
        right = Arrays.asList(new CigarElement(424, CigarOperator.H));
        data.add(new Object[]{alignment, left, middle, right});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "forCigarExtraction", groups = "sv")
    public void testExtractCigar(final AlignmentInterval interval, final List<CigarElement> expectedLeft,
                                 final List<CigarElement> expectedMiddle, final List<CigarElement> expectedRight) {

        final Tuple3<List<CigarElement>, List<CigarElement>, List<CigarElement>> x = SimpleStrandSwitchVariantDetector.splitCigarByLeftAndRightClipping(interval);
        Assert.assertEquals(x._1(), expectedLeft);
        Assert.assertEquals(x._2(), expectedMiddle);
        Assert.assertEquals(x._3(), expectedRight);
    }
}
