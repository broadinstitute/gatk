package org.broadinstitute.hellbender.utils.read;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;


public class FlowBasedReadUtilsUnitTest extends GATKBaseTest{
    @Test
    void testReadGroupParsing(){
        final String    testResourceDir = publicTestDir + "org/broadinstitute/hellbender/utils/read/flow/reads/";
        final String inputDir = testResourceDir + "/input/";

        final Path inputFile = FileSystems.getDefault().getPath(inputDir, "sample_mc.bam");
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(inputFile.toString()));
        SAMFileHeader header = reader.getFileHeader();
        SAMReadGroupRecord rg1 = header.getReadGroup("UGAv3-72");
        FlowBasedReadUtils.ReadGroupInfo frg1 = new FlowBasedReadUtils.ReadGroupInfo(rg1);
        assert(frg1.maxClass==12);
        assert(frg1.flowOrder.startsWith("TGCA"));
        assert(frg1.flowOrder.startsWith("TGCA"));
        assert(frg1.isFlowPlatform);

        SAMReadGroupRecord rg2 = header.getReadGroup("UGAv3-73");
        FlowBasedReadUtils.ReadGroupInfo frg2 = new FlowBasedReadUtils.ReadGroupInfo(rg2);
        assert(frg2.maxClass==20);
        assert(frg2.flowOrder.startsWith("TGCA"));
        assert(frg2.isFlowPlatform);
        SAMReadGroupRecord rg3 = header.getReadGroup("UGAv3-74");
        FlowBasedReadUtils.ReadGroupInfo frg3 = new FlowBasedReadUtils.ReadGroupInfo(rg3);
        assert(!frg3.isFlowPlatform);

    }

    @Test
    void testNonFlowReadGroupParsing(){
        final String    testResourceDir = publicTestDir + "org/broadinstitute/hellbender/utils/read/flow/reads/";
        final String inputDir = testResourceDir + "/input/";

        final Path inputFile = FileSystems.getDefault().getPath(inputDir, "non_flow_reads_with_tp.bam");
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(inputFile.toString()));
        SAMFileHeader header = reader.getFileHeader();
        SAMReadGroupRecord rg1 = header.getReadGroup("test");
        SAMRecord read = reader.iterator().next();
        GATKRead gread = new SAMRecordToGATKReadAdapter(read);
        assert(FlowBasedReadUtils.hasFlowTags(gread));
        FlowBasedReadUtils.ReadGroupInfo frg1 = FlowBasedReadUtils.getReadGroupInfo(header, gread);
        assert(!frg1.isFlowPlatform);
    }
    private GATKRead makeRead(final byte[] bases, final boolean isReverse) {

        byte[] quals = new byte[bases.length];

        final String cigar = String.format("%dM", bases.length);
        GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, cigar);
        read.setPosition("chr1",100);
        read.setIsReverseStrand(isReverse);
        return read;
    }

    @Test (dataProvider = "makeReads")
    public void testUncertainFlowTrimming(final GATKRead read, int nTrim, String uncertainFlowBase, final byte[] output, final int start, final int end) {
        FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();
        fbargs.flowNumUncertainFlows = nTrim;
        fbargs.flowFirstUncertainFlowBase = uncertainFlowBase;
        GATKRead trimmedRead = FlowBasedRead.hardClipUncertainBases(read, "TAGC", fbargs);
        Assert.assertEquals(trimmedRead.getBases(), output);
        Assert.assertEquals(trimmedRead.getStart(), start);
        Assert.assertEquals(trimmedRead.getEnd(), end);
    }


    @DataProvider(name="makeReads")
    public Object[][] makeMultipleReads(){
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, false), 4, "T", new byte[] {'G','A'}, 104, 105});
        tests.add(new Object[]{makeRead(new byte[]{'A','C','C','G','A','T'}, false), 4, "T", new byte[] {'G','A','T'}, 103, 105});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, true), 4, "T", new byte[] {'T','A','G','C'},100,103});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, false), 1, "T", new byte[] {'A','G','C','G','A'}, 101, 105});
        tests.add(new Object[]{makeRead(new byte[]{'A','C','C','G','A','T'}, false), 2, "T", new byte[] {'C','C','G','A','T'}, 101, 105});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, true), 10, "T", new byte[] {},0,0});

        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, false), 4, "A", new byte[] {'A','G','C','G','A'}, 101, 105});
        tests.add(new Object[]{makeRead(new byte[]{'A','C','C','G','A','T'}, false), 4, "A", new byte[] {'G','A','T'}, 103, 105});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, true), 4, "A", new byte[] {'T','A','G','C','G'},100,104});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, false), 1, "A", new byte[] {'T','A','G','C','G','A'}, 100, 105});
        tests.add(new Object[]{makeRead(new byte[]{'A','C','C','G','A','T'}, false), 2, "A", new byte[] {'C','C','G','A','T'}, 101, 105});
        tests.add(new Object[]{makeRead(new byte[]{'T','A','G','C','G','A'}, true), 10, "A", new byte[] {'T','A','G'},100,102});

        return tests.toArray(new Object[][]{});
    }


    @Test (dataProvider = "makeReads1")
    void testIsFlowPlatform(final SAMFileHeader hdr, final GATKRead read, final boolean answer){
        assert(FlowBasedReadUtils.isFlowPlatform(hdr, read) == answer);
    }

    @DataProvider(name="makeReads1")
    Object[][] makeReads1() {
        List<Object[]> tests = new ArrayList<>();
        SAMFileHeader hdr = new SAMFileHeader();

        GATKRead read1 = makeRead(new byte[]{'T','A','G','C','G','A'}, false);
        SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        rg1.setPlatform("ILLUMINA");
        read1.setReadGroup("rg1");
        hdr.addReadGroup(rg1);
        tests.add(new Object[]{hdr, read1, false});

        GATKRead read2 = makeRead(new byte[]{'T','A','G','C','G','A'}, false);
        read2.setAttribute("tp","blabla");
        SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");
        rg2.setPlatform("ILLUMINA");
        read2.setReadGroup("rg2");
        hdr.addReadGroup(rg2);
        tests.add(new Object[]{hdr, read2, false});

        GATKRead read3 = makeRead(new byte[]{'T','A','G','C','G','A'}, false);
        read3.setAttribute("tp",new byte[6]);
        SAMReadGroupRecord rg3 = new SAMReadGroupRecord("rg3");
        rg3.setPlatform("ULTIMA");
        read3.setReadGroup("rg3");
        hdr.addReadGroup(rg3);
        tests.add(new Object[]{hdr, read3, false});

        GATKRead read4 = makeRead(new byte[]{'T','A','G','C','G','A'}, false);
        read4.setAttribute("tp",new byte[6]);
        SAMReadGroupRecord rg4 = new SAMReadGroupRecord("rg4");
        rg4.setPlatform("ULTIMA");
        rg4.setAttribute("FO","TGCATGCA");
        read4.setReadGroup("rg4");
        hdr.addReadGroup(rg4);
        tests.add(new Object[]{hdr, read4, true});



        return tests.toArray(new Object[][]{});
    }
}
