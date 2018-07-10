package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.*;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class AssemblyBasedCallerUtilsUnitTest extends GATKBaseTest {
    private SAMFileHeader header;

    @BeforeClass
    public void init() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 100000000);
        SAMReadGroupRecord rg = new SAMReadGroupRecord("tumor");
        rg.setSample("tumor");
        header.addReadGroup(rg);
    }

    @Test
    public void testfinalizeRegion() {
        Assert.assertEquals(header.getSequenceIndex("1"), 0);
        final AssemblyRegion activeRegion = new AssemblyRegion(new SimpleInterval("1",42596728,42598843), 100, header);
        Assert.assertTrue(activeRegion.isActive());
        Assert.assertFalse(activeRegion.isFinalized());

        SAMLineParser parser = new SAMLineParser(header);
        List<GATKRead> reads = new LinkedList<GATKRead>();
        SAMRecord orgRead0 = parser.parseLine("HWI-ST807:461:C2P0JACXX:4:2204:18080:5857\t83\t1\t42596803\t39\t1S95M5S\t=\t42596891\t-7\tGAATCATCATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATAGAAT\tDAAAEDCFCCGEEDDBEDDDGCCDEDECDDFDCEECCFEECDCEDBCDBDBCC>DCECC>DBCDDBCBDDBCDDEBCCECC>DBCDBDBGC?FCCBDB>>?\tRG:Z:tumor");
        SAMRecord orgRead1 = parser.parseLine("HWI-ST807:461:C2P0JACXX:4:2204:18080:5857\t163\t1\t42596891\t39\t101M\t=\t42596803\t7\tCTCGAATGGAATCATTTTCTACTGGAAAGGAATGGAATCATCGCATAGAATCGAATGGAATTAACATGGAATGGAATCGAATGTAATCATCATCAAATGGA\t>@>:ABCDECCCEDCBBBDDBDDEBCCBEBBCBEBCBCDDCD>DECBGCDCF>CCCFCDDCBABDEDFCDCDFFDDDG?DDEGDDFDHFEGDDGECB@BAA\tRG:Z:tumor");

        reads.add(new SAMRecordToGATKReadAdapter(orgRead0.deepCopy()));
        reads.add(new SAMRecordToGATKReadAdapter(orgRead1.deepCopy()));
        activeRegion.addAll(reads);

        SampleList sampleList = SampleList.singletonSampleList("tumor");
        Byte minbq = 9;
        AssemblyBasedCallerUtils.finalizeRegion(activeRegion, false, false, minbq, header, sampleList);

        Assert.assertTrue(reads.get(0).convertToSAMRecord(header).equals(orgRead0));
        Assert.assertTrue(reads.get(1).convertToSAMRecord(header).equals(orgRead1));
    }
}
