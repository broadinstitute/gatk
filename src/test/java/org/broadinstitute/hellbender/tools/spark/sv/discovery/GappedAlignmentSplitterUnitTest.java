package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;
import java.util.stream.Collectors;

public class GappedAlignmentSplitterUnitTest {

    @Test(groups = "sv")
    public void testCompactifyNeighboringSoftClippings() {
        Assert.assertEquals(new Cigar(GappedAlignmentSplitter.compactifyNeighboringSoftClippings(TextCigarCodec.decode("1H2S3S4M5D6M7I8M9S10S11H").getCigarElements())),
                TextCigarCodec.decode("1H5S4M5D6M7I8M19S11H"));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_OneInsertion() {

        final Cigar cigar = TextCigarCodec.decode("56S27M15I32M21S");
        final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 158), 57, 130, cigar, true, 60, 0);

        final List<AlignedAssembly.AlignmentInterval> generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength())).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 126), 57, 83, TextCigarCodec.decode("56S27M68S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 127, 158), 99, 130, TextCigarCodec.decode("98S32M21S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));

    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_OneDeletion() {
        final Cigar cigar = TextCigarCodec.decode("2S205M2D269M77S");
        final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 575), 208, 476, cigar, true, 60, 0);

        final List<AlignedAssembly.AlignmentInterval> generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength())).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 304), 3, 207, TextCigarCodec.decode("2S205M346S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 307, 575), 208, 476, TextCigarCodec.decode("207S269M77S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_Complex() {

        final Cigar cigar = TextCigarCodec.decode("397S118M2D26M6I50M7I26M1I8M13D72M398S");
        final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 414), 398, 711, cigar, true, 60, 65);

        final List<AlignedAssembly.AlignmentInterval> generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength())).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 6);

        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 217), 398, 515, TextCigarCodec.decode("397S118M594S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 220, 245), 516, 541, TextCigarCodec.decode("515S26M568S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(2), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 246, 295), 548, 597, TextCigarCodec.decode("547S50M512S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(3), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 296, 321), 605, 630, TextCigarCodec.decode("604S26M479S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(4), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 322, 329), 632, 639, TextCigarCodec.decode("631S8M470S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(5), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 343, 414), 640, 711, TextCigarCodec.decode("639S72M398S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_GapSizeSensitivity() {

        final Cigar cigar = TextCigarCodec.decode("10M10D10M60I10M10I10M50D10M");
        final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 209), 1, 120, cigar, true, 60, 0);

        final List<AlignedAssembly.AlignmentInterval> generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, cigar.getReadLength())).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 3);
        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 129), 1, 20, TextCigarCodec.decode("10M10D10M100S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 130, 149), 81, 110, TextCigarCodec.decode("80S10M10I10M10S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(2), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 200, 209), 111, 120, TextCigarCodec.decode("110S10M"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_HardAndSoftClip() {

        final Cigar cigar = TextCigarCodec.decode("1H2S3M5I10M20D6M7S8H");
        final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 138), 4, 27, cigar, true, 60, 0);

        final List<AlignedAssembly.AlignmentInterval> generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength()+1+8)).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 100, 102), 4, 6, TextCigarCodec.decode("1H2S3M28S8H"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 103, 112), 12, 21, TextCigarCodec.decode("1H10S10M13S8H"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(2), new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 133, 138), 22, 27, TextCigarCodec.decode("1H20S6M7S8H"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_TerminalInsertionOperatorToSoftClip() {

        // beginning with 'I'
        Cigar cigar = TextCigarCodec.decode("10I10M5I10M");
        AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(new SimpleInterval("1", 101, 120), 11, 35, cigar, true, 60, 0);

        List<AlignedAssembly.AlignmentInterval> generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength())).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 101, 110), 11, 20, TextCigarCodec.decode("10S10M15S"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 111, 120), 26, 35, TextCigarCodec.decode("25S10M"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM));

        // ending with 'I'
        cigar = TextCigarCodec.decode("10M5I10M10I");
        alignmentInterval = new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 101, 120), 1, 25, cigar,
                true, 60, 0);

        generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength())).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 101, 110), 1, 10, TextCigarCodec.decode("10M25S"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 111, 120), 16, 25, TextCigarCodec.decode("15S10M10S"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_TerminalInsertionNeighboringClippings(){

        Cigar cigar = TextCigarCodec.decode("10H20S30I40M50I60S70H");
        AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 101, 140), 61, 100, cigar,
                true, 60, 0);
        List<AlignedAssembly.AlignmentInterval> generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength()+10+70)).collect(Collectors.toList());
        // no internal gap, so nothing should change
        Assert.assertEquals(generatedARList.size(), 1);
        Assert.assertEquals(generatedARList.get(0), alignmentInterval);

        cigar = TextCigarCodec.decode("10H20S30I40M5D15M50I60S70H");
        alignmentInterval = new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 101, 160), 61, 115, cigar,
                true, 60, 0);
        generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength()+10+70)).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);

        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 101, 140), 61, 100, TextCigarCodec.decode("10H50S40M125S70H"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 146, 160), 101, 115, TextCigarCodec.decode("10H90S15M110S70H"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_NegativeStrand() {
        // read data with AlignedAssembly.AlignmentInterval.toString():
        // 19149	contig-8	chrUn_JTFH01000557v1_decoy	21	1459	-	10S1044M122I395M75I	60	11	1646	200
        final Cigar cigar = TextCigarCodec.decode("10S1044M122I395M75I");
        final AlignedAssembly.AlignmentInterval alignmentInterval = new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 21, 1459), 11, 1646, cigar,
                false, 60, 200);
        final List<AlignedAssembly.AlignmentInterval> generatedARList = Utils.stream(GappedAlignmentSplitter.split(alignmentInterval, 1, cigar.getReadLength())).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.get(0), new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 416, 1459), 11, 1054, TextCigarCodec.decode("10S1044M592S"),
                false, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM));
        Assert.assertEquals(generatedARList.get(1), new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 21, 415), 1177, 1571, TextCigarCodec.decode("1176S395M75S"),
                false, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_ExpectException() {
        int exceptionCount = 0;

        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10S10D10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10D10S"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10D10H"));} catch (final Exception ex) {++exceptionCount;}

        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10S"));} catch (final Exception ex) {++exceptionCount;}

        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10M10D10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10M10I10S"));} catch (final Exception ex) {++exceptionCount;}

        // these 4 are fine now
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10S10I10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10I10S"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10I10H"));} catch (final Exception ex) {++exceptionCount;}

        // last two are valid
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10M10I10M10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10M10D10M10S"));} catch (final Exception ex) {++exceptionCount;}

        Assert.assertEquals(exceptionCount, 8);
    }

    private static Iterable<AlignedAssembly.AlignmentInterval> willThrowOnInvalidCigar(final Cigar cigar) throws GATKException {
        final AlignedAssembly.AlignmentInterval detailsDoesnotMatter = new AlignedAssembly.AlignmentInterval(
                new SimpleInterval("1", 1, 110), 21, 30,
                cigar,
                true, 60, 0);
        return GappedAlignmentSplitter.split(detailsDoesnotMatter, 1, cigar.getReadLength() + SVVariantDiscoveryUtils.getTotalHardClipping(cigar));
    }
}
