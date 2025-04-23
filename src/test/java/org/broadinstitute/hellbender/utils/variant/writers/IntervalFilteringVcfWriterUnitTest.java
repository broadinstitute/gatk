package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class IntervalFilteringVcfWriterUnitTest extends GATKBaseTest {


    @DataProvider
    public Object[][] getIntervalsAndMode(){
        final VariantContext noOverlap = new VariantContextBuilder("test", "1", 200, 300, Arrays.asList(Allele.create(Utils.repeatChars('A', 101), true), Allele.ALT_A)).make();
        final VariantContext contained =  new VariantContextBuilder("test", "1", 101, 104, Arrays.asList(Allele.create(Utils.repeatChars('A', 4), true), Allele.ALT_A)).make();
        final VariantContext overlaps =  new VariantContextBuilder("test", "1", 90, 120, Arrays.asList(Allele.create(Utils.repeatChars('A', 31), true), Allele.ALT_A)).make();
        final VariantContext startsIn =    new VariantContextBuilder("test", "1", 103, 140, Arrays.asList(Allele.create(Utils.repeatChars('A', 38), true), Allele.ALT_A)).make();
        final VariantContext endsIn =  new VariantContextBuilder("test", "1", 90, 103, Arrays.asList(Allele.create(Utils.repeatChars('A', 14), true), Allele.ALT_A)).make();
        final VariantContext anotherContig = new VariantContextBuilder("test", "2", 90, 140, Arrays.asList(Allele.create(Utils.repeatChars('A', 51), true), Allele.ALT_A)).make();
        final List<VariantContext> vcs = Arrays.asList(noOverlap, contained, overlaps, startsIn, endsIn, anotherContig);

        final SimpleInterval interval = new SimpleInterval("1", 100, 105);

        return new Object[][]{
                //                                                                      no overlap,   contained,   overlaps,   starts in,  ends in,  another contig
                {interval, vcs, IntervalFilteringVcfWriter.Mode.ANYWHERE,  new boolean[]{ true,         true,       true,       true,       true,     true}},
                {interval, vcs, IntervalFilteringVcfWriter.Mode.CONTAINED, new boolean[]{ false,        true,       false,      false,      false,    false}},
                {interval, vcs, IntervalFilteringVcfWriter.Mode.OVERLAPS,  new boolean[]{ false,        true,       true,       true,       true,     false}},
                {interval, vcs, IntervalFilteringVcfWriter.Mode.STARTS_IN, new boolean[]{ false,        true,       false,      true,       false,    false}},
                {interval, vcs, IntervalFilteringVcfWriter.Mode.ENDS_IN,   new boolean[]{ false,        true,       false,      false,      true,     false}},
        };
    }

    @Test(dataProvider = "getIntervalsAndMode")
    public void testModes(SimpleInterval interval, List<VariantContext> vcs, IntervalFilteringVcfWriter.Mode mode, boolean[] expected) {
        final OverlapDetector<SimpleInterval> detector = OverlapDetector.create(Collections.singletonList(interval));
        for(int i = 0; i < expected.length; i++){
            Assert.assertEquals(mode.test(detector,vcs.get(i)), expected[i], "mode " + mode + " mismatches at " + i);
        }
    }

    @Test
    public void testHeaderWriting() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final List<SimpleInterval> intervals = Arrays.asList(new SimpleInterval("1", 10, 100), new SimpleInterval("2", 100, 500));
        final IntervalFilteringVcfWriter writer = new IntervalFilteringVcfWriter(mockWriter, intervals, IntervalFilteringVcfWriter.Mode.OVERLAPS);
        writer.writeHeader(new VCFHeader());
        Assert.assertTrue(mockWriter.headerSet);
        Assert.assertTrue(mockWriter.headerWritten);
    }

    @Test
    public void testHeaderSetting(){
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final List<SimpleInterval> intervals = Arrays.asList(new SimpleInterval("1", 10, 100), new SimpleInterval("2", 100, 500));
        final IntervalFilteringVcfWriter writer = new IntervalFilteringVcfWriter(mockWriter, intervals, IntervalFilteringVcfWriter.Mode.OVERLAPS);
        writer.setHeader(new VCFHeader());
        Assert.assertTrue(mockWriter.headerSet);
        Assert.assertFalse(mockWriter.headerWritten);
    }

    @Test
    public void testClose() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final List<SimpleInterval> intervals = Arrays.asList(new SimpleInterval("1", 10, 100), new SimpleInterval("2", 100, 500));
        final IntervalFilteringVcfWriter writer = new IntervalFilteringVcfWriter(mockWriter, intervals, IntervalFilteringVcfWriter.Mode.OVERLAPS);
        writer.close();
        Assert.assertTrue(mockWriter.closed);
    }

    @Test
    public void testCheckError(){
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final List<SimpleInterval> intervals = Arrays.asList(new SimpleInterval("1", 10, 100), new SimpleInterval("2", 100, 500));
        final IntervalFilteringVcfWriter writer = new IntervalFilteringVcfWriter(mockWriter, intervals, IntervalFilteringVcfWriter.Mode.OVERLAPS);
        Assert.assertFalse(writer.checkError());
        mockWriter.error = true;
        Assert.assertTrue(writer.checkError());
    }

}