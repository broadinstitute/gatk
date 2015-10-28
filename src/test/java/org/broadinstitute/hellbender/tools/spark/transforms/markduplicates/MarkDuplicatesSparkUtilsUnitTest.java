package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.PairedEnds;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class MarkDuplicatesSparkUtilsUnitTest extends BaseTest {
    /*
        @DataProvider(name = "md")
        public Object[][] loadReads() {
            String dir = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/sam/MarkDuplicates/").getAbsolutePath();
            return new Object[][]{
                    {dir + "/example.chr1.1-1K.unmarkedDups.noDups.bam", 20, 0},
                    {dir + "/example.chr1.1-1K.unmarkedDups.bam", 90, 6},
                    {dir + "/example.chr1.1-1K.markedDups.bam", 90, 6},
            };
        }
        */
    @Test
    void handleFragmentsTest() {
        List<PairedEnds> pairedEnds = Collections.emptyList();
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        List<GATKRead> reads = MarkDuplicatesSparkUtils.handleFragments(pairedEnds, header);
        Assert.assertEquals(reads, Collections.emptyList());

        GATKRead myFakeRead = ArtificialReadUtils.createArtificialRead(header, "myFakeRead", 1, 200, 10);
        myFakeRead.setBaseQualities(Utils.dupBytes((byte) 30, 10));

        List<PairedEnds> singlePairedEnd = Collections.singletonList(PairedEnds.of(myFakeRead));
        List<GATKRead> noDupes = MarkDuplicatesSparkUtils.handleFragments(singlePairedEnd, header);
        Assert.assertEquals(noDupes, Collections.singletonList(myFakeRead));
    }

}