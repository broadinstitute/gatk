package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


public class SVSegmentUnitTest extends GATKBaseTest {

    @DataProvider(name="svSegmentPairs")
    public Object[][] getSVSegmentPairs() {
        return new Object[][] {
                // same SV type, different interval
                {
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 50, 150)),
                        1
                },
                // different SV type, use intervals
                {
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 50, 150)),
                        1
                },
                // same interval, use SV type
                {
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 100, 200)),
                        -1
                },
                // same segment
                {
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        0
                },
                // different contig, same coordinates/SV type
                {
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr2", 100, 200)),
                        -1
                },
                // different end coordinate only
                {
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 300)),
                        -1
                },
                // different start coordinate only
                {
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 99, 200)),
                        1
                },
        };
    }


    @Test(dataProvider = "svSegmentPairs")
    public void testCompareSVSegments(
            final SVSegment first,
            final SVSegment second,
            final int expectedComparisonValue)
    {
        final SAMSequenceDictionary dictionary =
                SVAnnotateUnitTest.createSequenceDictionary(Arrays.asList("chr1", "chr2", "chr3", "chr4"));
        final int actualComparisonValue = SVSegment.compareSVSegments(first, second, dictionary);
        Assert.assertEquals(actualComparisonValue, expectedComparisonValue);
    }


    @DataProvider(name="svSegmentList")
    public Object[][] getSVSegmentList() {
        return new Object[][] {
                {
                    Arrays.asList(
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 50, 150)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 50, 150)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr2", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 300)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 99, 200))

                    ), Arrays.asList(
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 50, 150)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 50, 150)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 99, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 200)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 100, 300)),
                        new SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr2", 100, 200))
                )
                },

        };
    }

    @Test(dataProvider = "svSegmentList")
    public void testSVSegmentComparator(
            final List<SVSegment> unsorted,
            final List<SVSegment> expectedOrder)
    {
        final SAMSequenceDictionary dictionary =
                SVAnnotateUnitTest.createSequenceDictionary(Arrays.asList("chr1", "chr2", "chr3", "chr4"));
        final List<SVSegment> actualOrder = unsorted.stream().sorted(SVSegment.getSVSegmentComparator(dictionary)).collect(Collectors.toList());
        Assert.assertEquals(actualOrder, expectedOrder);
    }
}
