package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

public class PrintDistantMatesUnitTest {
    @Test
    public void testPrintDistantMatesBasics() {
        final SAMFileHeader hdr =
                ArtificialReadUtils.createArtificialSamHeader(1, 1, 5000);
        final List<GATKRead> readPair =
                ArtificialReadUtils.createPair(hdr, "r1", 100, 0,
                1501, 3501, true, false);
        final GATKRead leftRead = readPair.get(0);
        final GATKRead rightRead = readPair.get(1);
        final GATKRead leftDistantMate = PrintDistantMates.doDistantMateAlterations(leftRead);
        Assert.assertEquals(leftDistantMate.getAssignedContig(), rightRead.getContig());
        Assert.assertEquals(leftDistantMate.getAssignedStart(), rightRead.getStart());
        Assert.assertTrue(PrintDistantMates.isDistantMate(leftDistantMate));
        Assert.assertEquals(leftRead, PrintDistantMates.undoDistantMateAlterations(leftDistantMate));
    }
}
