package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;


public class AssemblyAlignmentParserUnitTest extends BaseTest{

    @Test
    public void testGappedAlignmentBreaker_OneInsertion() {

        final Cigar cigar = TextCigarCodec.decode("56S27M15I32M21S");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 158), cigar, true, 60, 0, 57, 130);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1).spliterator(), false).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 126), TextCigarCodec.decode("56S27M68S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 57, 83));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 127, 158), TextCigarCodec.decode("98S32M21S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 99, 130));
    }

    @Test
    public void testGappedAlignmentBreaker_OneDeletion() {
        final Cigar cigar = TextCigarCodec.decode("2S205M2D269M77S");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 575), cigar, true, 60, 0, 208, 476);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1).spliterator(), false).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 304), TextCigarCodec.decode("2S205M346S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 3, 207));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 307, 575), TextCigarCodec.decode("207S269M77S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 208, 476));
    }

    @Test
    public void testGappedAlignmentBreaker_Complex() {

        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 414), TextCigarCodec.decode("397S118M2D26M6I50M7I26M1I8M13D72M398S"), true, 60, 65, 398, 711);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1).spliterator(), false).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 6);

        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 217), TextCigarCodec.decode("397S118M594S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 398, 515));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 220, 245), TextCigarCodec.decode("515S26M568S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 516, 541));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 246, 295), TextCigarCodec.decode("547S50M512S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 548, 597));
        Assert.assertEquals(generatedARList.get(3), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 296, 321), TextCigarCodec.decode("604S26M479S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 605, 630));
        Assert.assertEquals(generatedARList.get(4), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 322, 329), TextCigarCodec.decode("631S8M470S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 632, 639));
        Assert.assertEquals(generatedARList.get(5), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 343, 414), TextCigarCodec.decode("639S72M398S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 640, 711));
    }

    @Test
    public void testGappedAlignmentBreaker_Sensitivity() {

        final Cigar cigar = TextCigarCodec.decode("10M10D10M60I10M10I10M50D10M");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 209), cigar, true, 60, 0, 1, 120);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, SVConstants.CallingStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY).spliterator(), false).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 3);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 129), TextCigarCodec.decode("10M10D10M100S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 1, 20));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 130, 149), TextCigarCodec.decode("80S10M10I10M10S"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 81, 110));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 200, 209), TextCigarCodec.decode("110S10M"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 111, 120));
    }

    @Test
    public void testGappedAlignmentBreaker_HardAndSoftClip() {

        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 138), TextCigarCodec.decode("1H2S3M5I10M20D6M7S8H"), true, 60, 0, 4, 27);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1).spliterator(), false).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 102), TextCigarCodec.decode("1H2S3M28S8H"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 4, 6));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 103, 112), TextCigarCodec.decode("1H10S10M13S8H"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 12, 21));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 133, 138), TextCigarCodec.decode("1H20S6M7S8H"), true, 60, SVConstants.CallingStepConstants.MISSING_NM, 22, 27));
    }
}
