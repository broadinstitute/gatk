package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.testng.Assert.*;

public class GeneExpressionEvaluationUnitTest extends GATKBaseTest {

    @Test
    public void testInGoodPair() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 10000);
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header,"theRead", 1, 200, 100);


        read1.setMatePosition(read1.getContig(), 300);
        read1.setMateIsReverseStrand(true);
        read1.setIsProperlyPaired(true);
        read1.setAttribute(SAMTag.MQ.toString(), 20);
        read1.setFragmentLength(300);
        read1.setAttribute(SAMTag.MC.toString(), "100M");
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1));

        //not properly paired
        read1.setIsProperlyPaired(false);
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1));
        read1.setIsProperlyPaired(true);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1));

        //chimeric
        read1.setMatePosition(header.getSequenceDictionary().getSequence(0).getSequenceName(), 300);
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1));
        read1.setMatePosition(read1.getContig(), 300);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1));

        //same strand
        read1.setMateIsReverseStrand(false);

        read1.setMateIsReverseStrand(true);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1));

        //mate mapping quality too low
        read1.setAttribute(SAMTag.MQ.toString(), -1);
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1));
        read1.setAttribute(SAMTag.MQ.toString(), 20);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1));

        //outward facing
        read1.setMatePosition(read1.getContig(), 50);
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1));
        read1.setMatePosition(read1.getContig(), 250);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1));

        //swap strands and still pass
        read1.setIsReverseStrand(true);
        read1.setMateIsReverseStrand(false);
        read1.setPosition(read1.getContig(), 250);
        read1.setMatePosition(read1.getContig(), 200);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1));
    }

    @DataProvider(name = "testGetAlignmentIntervalsDataProvider")
    public Object[][] testGetAlignmentIntervalsDataProvider() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);

        final List<Object[]> examples = new ArrayList<>();
        final GATKRead splicedRead = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 100, 150);
        final String contig = splicedRead.getContig();
        splicedRead.setCigar("10M35N64M2I27M22N35M5D12M");
        splicedRead.setMatePosition(contig, 400);
        splicedRead.setMateIsReverseStrand(true);
        splicedRead.setIsProperlyPaired(true);
        splicedRead.setAttribute(SAMTag.MQ.toString(), 20);
        splicedRead.setAttribute(SAMTag.MC.toString(), "50M25N100M");
        examples.add(new Object[]{splicedRead, true, header, Arrays.asList(new Interval(contig, 100, 109),
                                                                           new Interval(contig, 145, 235),
                                                                           new Interval(contig, 258, 292),
                                                                           new Interval(contig, 298, 309),
                                                                           new Interval(contig, 400, 449),
                                                                           new Interval(contig, 475, 574)
                                                                          )
                                  }
                      );

        final GATKRead splicedReadImproperlyPaired = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 100, 150);
        splicedReadImproperlyPaired.setCigar("10M35N64M2I27M22N35M5D12M");
        splicedReadImproperlyPaired.setMatePosition(contig, 400);
        splicedReadImproperlyPaired.setMateIsReverseStrand(true);
        splicedReadImproperlyPaired.setIsProperlyPaired(false);
        splicedReadImproperlyPaired.setAttribute(SAMTag.MQ.toString(), 20);
        splicedReadImproperlyPaired.setAttribute(SAMTag.MC.toString(), "50M25N100M");
        examples.add(new Object[]{splicedReadImproperlyPaired, true, header, Arrays.asList(new Interval(contig, 100, 109),
                                                                           new Interval(contig, 145, 235),
                                                                           new Interval(contig, 258, 292),
                                                                           new Interval(contig, 298, 309)
                                                                           )
                                 }
        );

        final GATKRead splicedReadOverlappingMate = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 100, 150);
        splicedReadOverlappingMate.setCigar("10M35N64M2I27M22N35M5D12M");
        splicedReadOverlappingMate.setMatePosition(contig, 300);
        splicedReadOverlappingMate.setMateIsReverseStrand(true);
        splicedReadOverlappingMate.setIsProperlyPaired(true);
        splicedReadOverlappingMate.setAttribute(SAMTag.MQ.toString(), 20);
        splicedReadOverlappingMate.setAttribute(SAMTag.MC.toString(), "50M25N100M");
        examples.add(new Object[]{splicedReadOverlappingMate, true, header, Arrays.asList(new Interval(contig, 100, 109),
                new Interval(contig, 145, 235),
                new Interval(contig, 258, 292),
                new Interval(contig, 298, 349),
                new Interval(contig, 375, 474)
                )
                }
        );

        final GATKRead unsplicedRead = ArtificialReadUtils.createArtificialRead(header, "unsplicedRead", 0, 100, 150);
        unsplicedRead.setCigar("150M");
        unsplicedRead.setMatePosition(contig, 300);
        unsplicedRead.setMateIsReverseStrand(true);
        unsplicedRead.setIsProperlyPaired(true);
        unsplicedRead.setAttribute(SAMTag.MQ.toString(), 20);
        unsplicedRead.setAttribute(SAMTag.MC.toString(), "150M");
        unsplicedRead.setFragmentLength(350);
        examples.add(new Object[]{unsplicedRead, false, header, Collections.singletonList(new Interval(contig, 100, 449))}
        );

        final GATKRead unsplicedReadImproperlyPaired = ArtificialReadUtils.createArtificialRead(header, "unsplicedRead", 0, 100, 150);
        unsplicedReadImproperlyPaired.setCigar("150M");
        unsplicedReadImproperlyPaired.setMatePosition(contig, 300);
        unsplicedReadImproperlyPaired.setMateIsReverseStrand(true);
        unsplicedReadImproperlyPaired.setIsProperlyPaired(false);
        unsplicedReadImproperlyPaired.setAttribute(SAMTag.MQ.toString(), 20);
        unsplicedReadImproperlyPaired.setAttribute(SAMTag.MC.toString(), "150M");
        unsplicedReadImproperlyPaired.setFragmentLength(350);
        examples.add(new Object[]{unsplicedReadImproperlyPaired, false, header, Collections.singletonList(new Interval(contig, 100, 249))}
        );


        return examples.toArray(new Object[0][]);
    }

    @Test(dataProvider = "testGetAlignmentIntervalsDataProvider")
    public void testGetAlignmentIntervals(final GATKRead read, final boolean spliced, final SAMFileHeader header, final List<Interval> expectedAlignmentIntervals) {
        Assert.assertEquals(GeneExpressionEvaluation.getAlignmentIntervals(read, spliced, header), expectedAlignmentIntervals);
    }
}