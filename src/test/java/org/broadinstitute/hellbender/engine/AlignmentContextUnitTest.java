package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Map;

public final class AlignmentContextUnitTest extends GATKBaseTest {

    @Test
    public void test1Sample2Readgroups() throws Exception {
        final SAMReadGroupRecord readGroupOne = new SAMReadGroupRecord("rg1");
        readGroupOne.setSample("sample1");
        final SAMReadGroupRecord readGroupTwo = new SAMReadGroupRecord("rg2");
        readGroupTwo.setSample("sample1");

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        header.addReadGroup(readGroupOne);
        header.addReadGroup(readGroupTwo);

        final Locatable loc = new SimpleInterval("chr1", 1, 1);
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header,"read1",0,1,10);
        read1.setReadGroup(readGroupOne.getId());
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header,"read2",0,1,10);
        read2.setReadGroup(readGroupTwo.getId());
        final GATKRead read3 = ArtificialReadUtils.createArtificialRead(header,"read3",0,1,10);
        read3.setReadGroup(readGroupOne.getId());
        final GATKRead read4 = ArtificialReadUtils.createArtificialRead(header,"read4",0,1,10);
        read4.setReadGroup(readGroupTwo.getId());
        final GATKRead read5 = ArtificialReadUtils.createArtificialRead(header,"read5",0,1,10);
        read5.setReadGroup(readGroupTwo.getId());
        final GATKRead read6 = ArtificialReadUtils.createArtificialRead(header,"read6",0,1,10);
        read6.setReadGroup(readGroupOne.getId());
        final GATKRead read7 = ArtificialReadUtils.createArtificialRead(header,"read7",0,1,10);
        read7.setReadGroup(readGroupOne.getId());

        final ReadPileup pileup = new ReadPileup(loc, Arrays.asList(read1, read2, read3, read4, read5, read6, read7), 1);

        final AlignmentContext ac = new AlignmentContext(loc, pileup);
        Assert.assertSame(ac.getBasePileup(), pileup);
        Assert.assertEquals(ac.getContig(), loc.getContig());
        Assert.assertEquals(ac.getEnd(), loc.getEnd());
        Assert.assertEquals(ac.getLocation(), loc);
        Assert.assertEquals(ac.getPosition(), loc.getStart());
        Assert.assertEquals(ac.getStart(), loc.getStart());
        Assert.assertEquals(ac.hasPileupBeenDownsampled(), false);
        Assert.assertEquals(ac.size(), pileup.size());

        Assert.assertSame(ac.stratify(AlignmentContext.ReadOrientation.COMPLETE), ac, "Complete");

        final AlignmentContext acFWD = ac.stratify(AlignmentContext.ReadOrientation.FORWARD);
        Assert.assertEquals(acFWD.getLocation(), loc, "Forward Loc");
        Assert.assertEquals((Iterable<?>) acFWD.getBasePileup(), (Iterable<?>)pileup, "Forward Pileup");

        final AlignmentContext acREV = ac.stratify(AlignmentContext.ReadOrientation.REVERSE);
        AlignmentContext emptyAC= new AlignmentContext(loc, new ReadPileup(loc));
        Assert.assertEquals(acREV.getLocation(), loc, "Reverse Loc");
        Assert.assertEquals((Iterable<?>)acREV.getBasePileup(), (Iterable<?>)emptyAC.getBasePileup(), "Reverse pileup");
        Assert.assertNotNull(ac.toString());

        final Map<String, AlignmentContext> bySample = ac.splitContextBySampleName(header);
        Assert.assertEquals(bySample.size(), 1);
        Assert.assertEquals(bySample.keySet(), ReadUtils.getSamplesFromHeader(header));
        final AlignmentContext firstAC = bySample.values().iterator().next();
        Assert.assertEquals(firstAC.getLocation(), ac.getLocation());
        Assert.assertEquals(firstAC.getBasePileup(), ac.getBasePileup());

        final Map<String, AlignmentContext> bySampleAssume1 = ac.splitContextBySampleName("sample1", header);
        Assert.assertEquals(bySampleAssume1.size(), 1);
        Assert.assertEquals(bySampleAssume1.keySet(), ReadUtils.getSamplesFromHeader(header));
        final AlignmentContext firstACAssume1 = bySampleAssume1.values().iterator().next();
        Assert.assertEquals(firstACAssume1.getLocation(), ac.getLocation());
        Assert.assertEquals(firstACAssume1.getBasePileup(), ac.getBasePileup());

        final Map<String, AlignmentContext> stringAlignmentContextMap = AlignmentContext.splitContextBySampleName(pileup, header);
        Assert.assertEquals(stringAlignmentContextMap.keySet(), Collections.singleton("sample1"));
        Assert.assertEquals(stringAlignmentContextMap.get("sample1").getLocation(), loc);
        Assert.assertEquals(stringAlignmentContextMap.get("sample1").getBasePileup(), pileup);
    }


    @Test(expectedExceptions = UserException.ReadMissingReadGroup.class)
    public void testNoSample() throws Exception {
        final SAMReadGroupRecord readGroupOne = new SAMReadGroupRecord("rg1");

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        header.addReadGroup(readGroupOne);

        final Locatable loc = new SimpleInterval("chr1", 1, 1);
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header,"read1",0,1,10);
        read1.setReadGroup(readGroupOne.getId());

        final ReadPileup pileup = new ReadPileup(loc, Arrays.asList(read1), 1);
        final AlignmentContext ac = new AlignmentContext(loc, pileup);
        ac.splitContextBySampleName(header);

    }
}
