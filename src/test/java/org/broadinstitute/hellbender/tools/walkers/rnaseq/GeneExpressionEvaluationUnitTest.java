package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3BaseData;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GeneExpressionEvaluationUnitTest extends GATKBaseTest {
    private static final double DEFAULT_LENIENCE = 0.00001;

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
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));

        //not properly paired
        read1.setIsProperlyPaired(false);
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));
        read1.setIsProperlyPaired(true);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));

        //chimeric
        read1.setMatePosition(header.getSequenceDictionary().getSequence(0).getSequenceName(), 300);
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));
        read1.setMatePosition(read1.getContig(), 300);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));

        //strand orientations
        //reads are +-
        //both on transcription strand or neither on transcription strand, bad
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD));
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE));
        //one on transcription strand, other not, good
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD));

        //reads are ++
        read1.setMateIsReverseStrand(false);
        //both on transcription strand or neither on transcription strand, good
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD));
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE));
        //one on transcription strand, other not, bad
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD));

        //reads are -+
        read1.setIsReverseStrand(true);
        //both on transcription strand or neither on transcription strand, bad
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD));
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE));
        //one on transcription strand, other not, good
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD));

        //reads are --
        read1.setMateIsReverseStrand(true);
        //both on transcription strand or neither on transcription strand, good
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD));
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE));
        //one on transcription strand, other not, bad
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD));

        //mate mapping quality too low
        read1.setIsReverseStrand(false);
        read1.setAttribute(SAMTag.MQ.toString(), -1);
        Assert.assertFalse(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));
        read1.setAttribute(SAMTag.MQ.toString(), 20);
        Assert.assertTrue(GeneExpressionEvaluation.inGoodPair(read1, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE));
    }

    @DataProvider(name = "testGetAlignmentIntervalsDataProvider")
    public Object[][] testGetAlignmentIntervalsDataProvider() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000);

        final List<Object[]> examples = new ArrayList<>();
        final GATKRead splicedRead = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 100, 150);
        final String contig = splicedRead.getContig();
        splicedRead.setCigar("10M35N64M2I27M22N35M5D12M");
        pairRead(splicedRead, contig, 400, "50M25N100M");
        final List<Interval> expectedIntervalsSplicedRead = Arrays.asList(new Interval(contig, 100, 109),
                new Interval(contig, 145, 235),
                new Interval(contig, 258, 292),
                new Interval(contig, 298, 309),
                new Interval(contig, 400, 449),
                new Interval(contig, 475, 574)
        );
        examples.add(new Object[]{splicedRead, true, expectedIntervalsSplicedRead});

        final GATKRead splicedReadReverseStrand = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 400, 150);
        splicedReadReverseStrand.setCigar("50M25N100M");
        splicedReadReverseStrand.setIsReverseStrand(true);
        pairRead(splicedReadReverseStrand, contig, 100, "10M35N64M2I27M22N35M5D12M");
        examples.add(new Object[]{splicedReadReverseStrand, true, expectedIntervalsSplicedRead});

        final GATKRead splicedReadImproperlyPaired = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 100, 150);
        splicedReadImproperlyPaired.setCigar("10M35N64M2I27M22N35M5D12M");
        pairRead(splicedReadImproperlyPaired, contig, 400, "50M25N100M");
        splicedReadImproperlyPaired.setIsProperlyPaired(false);
        examples.add(new Object[]{splicedReadImproperlyPaired, true, Arrays.asList(new Interval(contig, 100, 109),
                                                                           new Interval(contig, 145, 235),
                                                                           new Interval(contig, 258, 292),
                                                                           new Interval(contig, 298, 309)
                                                                           )
                                 }
        );

        final GATKRead splicedReadImproperlyPairedReverseStrand = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 400, 150);
        splicedReadImproperlyPairedReverseStrand.setCigar("50M25N100M");
        splicedReadImproperlyPairedReverseStrand.setIsReverseStrand(true);
        pairRead(splicedReadImproperlyPairedReverseStrand, contig, 100, "10M35N64M2I27M22N35M5D12M");
        splicedReadImproperlyPairedReverseStrand.setIsProperlyPaired(false);
        examples.add(new Object[]{splicedReadImproperlyPairedReverseStrand, true, Arrays.asList(new Interval(contig, 400, 449), new Interval(contig, 475, 574))}
        );

        final GATKRead splicedReadOverlappingMate = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 100, 150);
        final List<Interval> expectedIntervalsSplicedReadOverlappingMate = Arrays.asList(new Interval(contig, 100, 109),
                new Interval(contig, 145, 235),
                new Interval(contig, 258, 292),
                new Interval(contig, 298, 349),
                new Interval(contig, 375, 474)
        );

        splicedReadOverlappingMate.setCigar("10M35N64M2I27M22N35M5D12M");
        pairRead(splicedReadOverlappingMate, contig, 300, "50M25N100M");
        examples.add(new Object[]{splicedReadOverlappingMate, true, expectedIntervalsSplicedReadOverlappingMate});

        final GATKRead splicedReadOverlappingMateReverseStrand = ArtificialReadUtils.createArtificialRead(header, "splicedRead", 0, 300, 150);
        splicedReadOverlappingMateReverseStrand.setCigar("50M25N100M");
        splicedReadOverlappingMateReverseStrand.setIsReverseStrand(true);
        pairRead(splicedReadOverlappingMateReverseStrand, contig, 100, "10M35N64M2I27M22N35M5D12M");
        examples.add(new Object[]{splicedReadOverlappingMateReverseStrand, true, expectedIntervalsSplicedReadOverlappingMate});

        final GATKRead unsplicedRead = ArtificialReadUtils.createArtificialRead(header, "unsplicedRead", 0, 100, 150);
        final List<Interval> expectedIntervalsUnsplicedRead = Collections.singletonList(new Interval(contig, 100, 449));
        unsplicedRead.setCigar("150M");
        pairRead(unsplicedRead, contig, 300, "150M");
        unsplicedRead.setFragmentLength(350);
        examples.add(new Object[]{unsplicedRead, false, expectedIntervalsUnsplicedRead}
        );

        final GATKRead unsplicedReadReverseStrand = ArtificialReadUtils.createArtificialRead(header, "unsplicedRead", 0, 300, 150);
        unsplicedReadReverseStrand.setCigar("150M");
        unsplicedReadReverseStrand.setIsReverseStrand(true);
        pairRead(unsplicedReadReverseStrand, contig, 100, "150M");
        unsplicedReadReverseStrand.setFragmentLength(-350);
        examples.add(new Object[]{unsplicedReadReverseStrand, false, expectedIntervalsUnsplicedRead}
        );

        final GATKRead unsplicedReadImproperlyPaired = ArtificialReadUtils.createArtificialRead(header, "unsplicedRead", 0, 100, 150);
        unsplicedReadImproperlyPaired.setCigar("150M");
        pairRead(unsplicedReadImproperlyPaired, contig, 300, "150M");
        unsplicedReadImproperlyPaired.setIsProperlyPaired(false);
        unsplicedReadImproperlyPaired.setFragmentLength(350);
        examples.add(new Object[]{unsplicedReadImproperlyPaired, false, Collections.singletonList(new Interval(contig, 100, 249))});

        final GATKRead unsplicedReadImproperlyPairedReverseStrand = ArtificialReadUtils.createArtificialRead(header, "unsplicedRead", 0, 300, 150);
        unsplicedReadImproperlyPairedReverseStrand.setCigar("150M");
        unsplicedReadImproperlyPairedReverseStrand.setIsReverseStrand(true);
        pairRead(unsplicedReadImproperlyPairedReverseStrand, contig, 100, "150M");
        unsplicedReadImproperlyPairedReverseStrand.setIsProperlyPaired(false);
        unsplicedReadImproperlyPairedReverseStrand.setFragmentLength(-350);
        examples.add(new Object[]{unsplicedReadImproperlyPairedReverseStrand, false, Collections.singletonList(new Interval(contig, 300, 449))});


        return examples.toArray(new Object[0][]);
    }

    private void pairRead(final GATKRead read, final String mateContig, final int mateStart, final String mateCigar) {
        read.setIsProperlyPaired(true);
        read.setMatePosition(mateContig, mateStart);
        read.setMateIsReverseStrand(!read.isReverseStrand());
        read.setAttribute(SAMTag.MQ.toString(), 20);
        read.setAttribute(SAMTag.MC.toString(), mateCigar);
    }

    @Test(dataProvider = "testGetAlignmentIntervalsDataProvider")
    public void testGetAlignmentIntervals(final GATKRead read, final boolean spliced, final List<Interval> expectedAlignmentIntervals) {
        Assert.assertEquals(GeneExpressionEvaluation.getAlignmentIntervals(read, !spliced, 0, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE), expectedAlignmentIntervals);
    }

    @DataProvider(name = "testSenseAntisenseDataProvider")
    public Object[][] testSenseAntisenseDataProvider() {
        final Gff3BaseData forwardFeature = new Gff3BaseData("1", ".", "gene", 1, 100, -1d, Strand.POSITIVE, 0, Collections.emptyMap());
        final Gff3BaseData reverseFeature = new Gff3BaseData("1", ".", "gene", 1, 100, -1d, Strand.NEGATIVE, 0, Collections.emptyMap());
        final Gff3BaseData unstrandedFeature = new Gff3BaseData("1", ".", "gene", 1, 100, -1d, Strand.NONE, 0, Collections.emptyMap());

        final List<Object[]> examples = new ArrayList<>();
        final GATKRead read1Forward = getReadForStrandTest(false, true);
        examples.add(new Object[]{read1Forward, forwardFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, true});
        examples.add(new Object[]{read1Forward, forwardFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, true});
        examples.add(new Object[]{read1Forward, forwardFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, false});
        examples.add(new Object[]{read1Forward, forwardFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, false});

        examples.add(new Object[]{read1Forward, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, true});
        examples.add(new Object[]{read1Forward, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, true});
        examples.add(new Object[]{read1Forward, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, false});
        examples.add(new Object[]{read1Forward, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, false});

        examples.add(new Object[]{read1Forward, reverseFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, false});
        examples.add(new Object[]{read1Forward, reverseFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, false});
        examples.add(new Object[]{read1Forward, reverseFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, true});
        examples.add(new Object[]{read1Forward, reverseFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, true});

        final GATKRead read1Reverse = getReadForStrandTest(true, true);
        examples.add(new Object[]{read1Reverse, forwardFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, false});
        examples.add(new Object[]{read1Reverse, forwardFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, false});
        examples.add(new Object[]{read1Reverse, forwardFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, true});
        examples.add(new Object[]{read1Reverse, forwardFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, true});

        examples.add(new Object[]{read1Reverse, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, false});
        examples.add(new Object[]{read1Reverse, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, false});
        examples.add(new Object[]{read1Reverse, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, true});
        examples.add(new Object[]{read1Reverse, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, true});

        examples.add(new Object[]{read1Reverse, reverseFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, true});
        examples.add(new Object[]{read1Reverse, reverseFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, true});
        examples.add(new Object[]{read1Reverse, reverseFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, false});
        examples.add(new Object[]{read1Reverse, reverseFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, false});

        final GATKRead read2Forward = getReadForStrandTest(false, false);
        examples.add(new Object[]{read2Forward, forwardFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, true});
        examples.add(new Object[]{read2Forward, forwardFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, false});
        examples.add(new Object[]{read2Forward, forwardFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, true});
        examples.add(new Object[]{read2Forward, forwardFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, false});

        examples.add(new Object[]{read2Forward, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, true});
        examples.add(new Object[]{read2Forward, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, false});
        examples.add(new Object[]{read2Forward, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, true});
        examples.add(new Object[]{read2Forward, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, false});

        examples.add(new Object[]{read2Forward, reverseFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, false});
        examples.add(new Object[]{read2Forward, reverseFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, true});
        examples.add(new Object[]{read2Forward, reverseFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, false});
        examples.add(new Object[]{read2Forward, reverseFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, true});

        final GATKRead read2Reverse = getReadForStrandTest(true, false);
        examples.add(new Object[]{read2Reverse, forwardFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, false});
        examples.add(new Object[]{read2Reverse, forwardFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, true});
        examples.add(new Object[]{read2Reverse, forwardFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, false});
        examples.add(new Object[]{read2Reverse, forwardFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, true});

        examples.add(new Object[]{read2Reverse, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, false});
        examples.add(new Object[]{read2Reverse, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, true});
        examples.add(new Object[]{read2Reverse, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, false});
        examples.add(new Object[]{read2Reverse, unstrandedFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, true});

        examples.add(new Object[]{read2Reverse, reverseFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_FORWARD, true});
        examples.add(new Object[]{read2Reverse, reverseFeature, GeneExpressionEvaluation.ReadStrands.FORWARD_REVERSE, false});
        examples.add(new Object[]{read2Reverse, reverseFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_FORWARD, true});
        examples.add(new Object[]{read2Reverse, reverseFeature, GeneExpressionEvaluation.ReadStrands.REVERSE_REVERSE, false});
        return examples.toArray(new Object[0][]);
    }

    private GATKRead getReadForStrandTest(final boolean isReverseStrand, final boolean isFirstOfPair) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead("150M");
        read.setIsReverseStrand(isReverseStrand);
        if (isFirstOfPair) {
            read.setIsFirstOfPair();
        } else {
            read.setIsSecondOfPair();
        }

        return read;
    }

    @Test(dataProvider = "testSenseAntisenseDataProvider")
    public void testSenseAntisense(final GATKRead read, final Gff3BaseData feature, final GeneExpressionEvaluation.ReadStrands readStrands, final boolean expectedSense) {
        final boolean actualSense = readStrands.isSense(read, feature);
        Assert.assertEquals(actualSense, expectedSense);
    }


    @DataProvider(name = "testMultiOverlapMethodDataProvider")
    public Object[][] testMultiOverlapMethodDataProvider() {
        /*
        first, build overlap detector.  three genes in play

                           150-250         350-450
         gene1             --------        --------
                               225-300                500-600
         gene2                 ---------              --------
                      100-200               400-410
         gene3        -------                 ---
        */

        final Gff3BaseData gene1 = new Gff3BaseData("theContig", ".", "gene", 150, 450, -1d, Strand.POSITIVE, 0, Collections.singletonMap("ID", Collections.singletonList("gene1")));
        final Gff3BaseData gene2 = new Gff3BaseData("theContig", ".", "gene", 225, 600, -1d, Strand.POSITIVE, 0, Collections.singletonMap("ID", Collections.singletonList("gene2")));
        final Gff3BaseData gene3 = new Gff3BaseData("theContig", ".", "gene", 100, 410, -1d, Strand.POSITIVE, 0, Collections.singletonMap("ID", Collections.singletonList("gene3")));

        final OverlapDetector<Pair<Gff3BaseData, Interval>> featureOverlapDetector = new OverlapDetector<>(0,0);
        //gene1 overlaps
        for (final Interval interval : Arrays.asList(new Interval("theContig", 150, 250), new Interval("theContig", 350, 450))) {
            featureOverlapDetector.addLhs(Pair.of(gene1, interval), interval);
        }

        //gene2 overlaps
        for (final Interval interval : Arrays.asList(new Interval("theContig", 225, 300), new Interval("theContig", 500, 600))) {
            featureOverlapDetector.addLhs(Pair.of(gene2, interval), interval);
        }

        //gene3 overlaps
        for (final Interval interval : Arrays.asList(new Interval("theContig", 100, 200), new Interval("theContig", 400, 410))) {
            featureOverlapDetector.addLhs(Pair.of(gene3, interval), interval);
        }

        final List<Object[]> examples = new ArrayList<>();

        /*
        alignments overlap none of genes, expect empty map

                           150-250         350-450
         gene1             --------        --------
                               225-300                500-600
         gene2                 ---------              --------
                      100-200               400-410
         gene3        -------                 ---
               20-75                   325-340                           700-750
    alignments  ***                      **                              ********
        */
        final List<Interval> noOverlapAlignments = Arrays.asList(new Interval("theContig", 20, 75), new Interval("theContig", 325, 340), new Interval("theContig", 700,750));
        examples.add(new Object[]{GeneExpressionEvaluation.MultiOverlapMethod.EQUAL, noOverlapAlignments, featureOverlapDetector, Collections.emptyMap()});
        examples.add(new Object[]{GeneExpressionEvaluation.MultiOverlapMethod.PROPORTIONAL, noOverlapAlignments, featureOverlapDetector, Collections.emptyMap()});

        /*
        alignments overlap only gene1, but not fully contained in gene1 exons

                           150-250         350-450
         gene1             --------        --------
                               225-300                500-600
         gene2                 ---------              --------
                      100-200               400-410
         gene3        -------                 ---
                                       325-375  440-460
    alignments                           ***      **
        */
        final List<Interval> overlapGene1Alignments = Arrays.asList(new Interval("theContig", 325, 375), new Interval("theContig", 440, 460));
        //for EQUAL method, gene1 has weight 1
        final Map<Gff3BaseData, Double> overlapGene1EqualWeights = Collections.singletonMap(gene1, (double)1);
        examples.add(new Object[]{GeneExpressionEvaluation.MultiOverlapMethod.EQUAL, overlapGene1Alignments, featureOverlapDetector, overlapGene1EqualWeights});

        //for PROPORTIONAL method, alignments are 72 bases, with 37 overlapping gene1.  So gene1 has weight 37/72
        final Map<Gff3BaseData, Double> overlapGene1ProportionalWeights = Collections.singletonMap(gene1, (double)37/(double)72);
        examples.add(new Object[]{GeneExpressionEvaluation.MultiOverlapMethod.PROPORTIONAL, overlapGene1Alignments, featureOverlapDetector, overlapGene1ProportionalWeights});

        /*
        alignments overlap gene1 and gene3

                           150-250         350-450
         gene1             --------        --------
                               225-300                500-600
         gene2                 ---------              --------
                      100-200               400-410
         gene3        -------                 ---
                       125-200              400-410
    alignments          *****                 ***
        */
        final List<Interval> overlapGene1AndGene3Alignments = Arrays.asList(new Interval("theContig", 125, 200), new Interval("theContig", 400, 410));
        //for EQUAL method, gene1 and gene 3 both have weight 1/2
        final Map<Gff3BaseData, Double> overlapGene1Gene3EqualWeights = new HashMap<>();
        Arrays.asList(gene1,gene3).forEach(g -> overlapGene1Gene3EqualWeights.put(g, (double)1/(double)2));


        examples.add(new Object[]{GeneExpressionEvaluation.MultiOverlapMethod.EQUAL, overlapGene1AndGene3Alignments, featureOverlapDetector, overlapGene1Gene3EqualWeights});
        /* for PROPORTIONAL METHOD
        gene1 alignment bases: 51 + 11 = 62
        gene2 alignment bases: 0
        gene3 alignment bases: 76 + 11 = 87
        off gene alignment bases: 0
        gene1 weight = 62/(62+87)
        gene3 weight = 87/(62+87)
         */
        final Map<Gff3BaseData, Double> overlapGene1Gene3ProportionalWeights = new HashMap<>();
        overlapGene1Gene3ProportionalWeights.put(gene1, (double)62/(double)(62 + 87));
        overlapGene1Gene3ProportionalWeights.put(gene3, (double)87/(double)(62 + 87));
        examples.add(new Object[]{GeneExpressionEvaluation.MultiOverlapMethod.PROPORTIONAL, overlapGene1AndGene3Alignments, featureOverlapDetector, overlapGene1Gene3ProportionalWeights});

        /*
        alignments overlap all three genes, and extend into region covered my none of genes

                           150-250         350-450
         gene1             --------        --------
                               225-300                500-600
         gene2                 ---------              --------
                      100-200               400-410
         gene3        -------                 ---
                           175-250              400-525
    alignments              *******             *******
        */

        final List<Interval> overlapAll3Alignments = Arrays.asList(new Interval("theContig", 175, 250), new Interval("theContig", 400, 525));
        //for EQUAL method, all three genes have wight 1/3
        final Map<Gff3BaseData, Double> overlapAll3EqualWeights = new HashMap<>();
        Arrays.asList(gene1, gene2, gene3).forEach(g -> overlapAll3EqualWeights.put(g, (double)1/(double)3));
        examples.add(new Object[]{GeneExpressionEvaluation.MultiOverlapMethod.EQUAL, overlapAll3Alignments, featureOverlapDetector, overlapAll3EqualWeights});
        /* for PROPORTIONAL METHOD
        gene1 alignment bases: 76 + 51 = 127
        gene2 alignment bases: 26 + 26 = 52
        gene3 alignment bases: 26 + 11 = 37
        off gene alignment bases: 49
        gene1 weight = 127/(127 + 52 + 37 + 49)
        gene2 weight = 52/(127 + 52 + 37 + 49)
        gene3 weight = 37/(127 + 52 + 37 + 49)
         */
        final Map<Gff3BaseData, Double> overlapAll3ProportionalWeights = new HashMap<>();
        overlapAll3ProportionalWeights.put(gene1, (double)127/(double)(127 + 52 + 37 + 49));
        overlapAll3ProportionalWeights.put(gene2, (double)52/(double)(127 + 52 + 37 + 49));
        overlapAll3ProportionalWeights.put(gene3, (double)37/(double)(127 + 52 + 37 + 49));
        examples.add(new Object[]{GeneExpressionEvaluation.MultiOverlapMethod.PROPORTIONAL, overlapAll3Alignments, featureOverlapDetector, overlapAll3ProportionalWeights});


        return examples.toArray(new Object[0][]);
    }

    @Test(dataProvider = "testMultiOverlapMethodDataProvider")
    public void testMultiOverlapMethod(final GeneExpressionEvaluation.MultiOverlapMethod multiOverlapMethod, final List<Interval> alignmentIntervals, final OverlapDetector<Pair<Gff3BaseData, Interval>> featureOverlapDetector,
                                       final Map<Gff3BaseData, Double> expectedWeights) {
        final Map<Gff3BaseData, Double> actualWeights = multiOverlapMethod.getWeights(alignmentIntervals, featureOverlapDetector);
        assertWeightMapsEquivalent(actualWeights, expectedWeights);
    }

    @DataProvider(name = "testMultiMapMethodDataProvider")
    public Object[][] testMultiMapMethodDataProvider() {
        final Gff3BaseData gene1 = new Gff3BaseData("theContig", ".", "gene", 150, 450, -1d, Strand.POSITIVE, 0, Collections.singletonMap("ID", Collections.singletonList("gene1")));
        final Gff3BaseData gene2 = new Gff3BaseData("theContig", ".", "gene", 225, 600, -1d, Strand.POSITIVE, 0, Collections.singletonMap("ID", Collections.singletonList("gene2")));
        final Gff3BaseData gene3 = new Gff3BaseData("theContig", ".", "gene", 100, 410, -1d, Strand.POSITIVE, 0, Collections.singletonMap("ID", Collections.singletonList("gene3")));

        final Map<Gff3BaseData, Double> previousWeights = new HashMap<>();
        previousWeights.put(gene1, (double)47/(double)3);
        previousWeights.put(gene2, (double)93/(double)7);
        previousWeights.put(gene3, (double)7653/(double)2);

        final List<Object[]> examples = new ArrayList<>();


        examples.add(new Object[]{GeneExpressionEvaluation.MultiMapMethod.IGNORE, 1, previousWeights, previousWeights});
        examples.add(new Object[]{GeneExpressionEvaluation.MultiMapMethod.IGNORE, 3, previousWeights, Collections.emptyMap()});
        examples.add(new Object[]{GeneExpressionEvaluation.MultiMapMethod.IGNORE, 7, previousWeights, Collections.emptyMap()});

        examples.add(new Object[]{GeneExpressionEvaluation.MultiMapMethod.EQUAL, 1, previousWeights, previousWeights});

        final Map<Gff3BaseData, Double> expectedWeightsEqualNhits3 = new HashMap<>();
        expectedWeightsEqualNhits3.put(gene1, (double)47/(double)(3*3));
        expectedWeightsEqualNhits3.put(gene2, (double)93/(double)(7*3));
        expectedWeightsEqualNhits3.put(gene3, (double)7653/(double)(2*3));
        examples.add(new Object[]{GeneExpressionEvaluation.MultiMapMethod.EQUAL, 3, previousWeights, expectedWeightsEqualNhits3});

        final Map<Gff3BaseData, Double> expectedWeightsEqualNhits7 = new HashMap<>();
        expectedWeightsEqualNhits7.put(gene1, (double)47/(double)(3*7));
        expectedWeightsEqualNhits7.put(gene2, (double)93/(double)(7*7));
        expectedWeightsEqualNhits7.put(gene3, (double)7653/(double)(2*7));
        examples.add(new Object[]{GeneExpressionEvaluation.MultiMapMethod.EQUAL, 7, previousWeights, expectedWeightsEqualNhits7});

        return examples.toArray(new Object[0][]);
    }

    @Test(dataProvider = "testMultiMapMethodDataProvider")
    public void testMultiMapMethod(final GeneExpressionEvaluation.MultiMapMethod multiMapMethod, final int nHits, final Map<Gff3BaseData, Double> previousWeights, final Map<Gff3BaseData, Double> expectedWeights) {
        final Map<Gff3BaseData, Double> actualWeights = multiMapMethod.getWeights(nHits, previousWeights);
        assertWeightMapsEquivalent(actualWeights, expectedWeights);
    }

    private <T> void assertWeightMapsEquivalent(final Map<T, Double> actualMap, final Map<T, Double> expectedMap) {
        Assert.assertEquals(actualMap.keySet(), expectedMap.keySet());
        for (final Map.Entry<T, Double> entry : actualMap.entrySet()) {
            final T key = entry.getKey();
            final double actualValue = entry.getValue();
            final double expectedValue = expectedMap.get(key);

            //should never have non-zero weights
            Assert.assertTrue(actualValue >= 0);
            Assert.assertTrue(expectedValue >= 0);

            if (actualValue + expectedValue > 0) {
                Assert.assertTrue(2 * Math.abs(actualValue - expectedValue) / (actualValue + expectedValue) <= DEFAULT_LENIENCE, " for key: " + key.toString() + " actual = " + actualValue + ", expected = " + expectedValue);
            }
        }
    }
}
