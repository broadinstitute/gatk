package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMTag;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

public class MoleculeIDUnitTest extends BaseTest {
    @Test
    public void testBasicOperations(){
        final int moleculeNumber = 131;
        final String strand = "A";
        final String samField = moleculeNumber + "/" + strand;
        final GATKRead read = ArtificialReadUtils.createSamBackedRead( "read", "3", 1000, 1151);
        read.setAttribute(SAMTag.MI.name(), samField);
        final MoleculeID id = new MoleculeID(read);
        Assert.assertEquals(id.getMoleculeNumber(), moleculeNumber);
        Assert.assertEquals(id.getStrand(), strand);
        Assert.assertEquals(id.getSAMField(), samField);
    }

    @Test
    public void testCountStrand(){
        final int numReadPairs = 50;
        final List<GATKRead> reads = ReadsWithSameUMIUnitTest.makeReads(numReadPairs, "3", 10, "read");
        Assert.assertEquals((int) MoleculeID.countStrands(reads).getLeft(), 50);
        Assert.assertEquals((int) MoleculeID.countStrands(reads).getRight(), 50);
    }
}