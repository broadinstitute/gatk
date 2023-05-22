package org.broadinstitute.hellbender.utils.read;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;


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
}
