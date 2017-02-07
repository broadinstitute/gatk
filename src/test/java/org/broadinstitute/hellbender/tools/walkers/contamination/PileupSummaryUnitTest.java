package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

/**
 * Created by David Benjamin on 2/15/17.
 */
public class PileupSummaryUnitTest {

    @Test
    public void test() throws IOException {
        final String contig = "chr1";
        final int position = 1;
        final int refCount = 20;
        final int altCount = 10;
        final int otherAltCount = 2;
        final double alleleFrequency = 0.3;

        final List<PileupSummary> ps = Arrays.asList(new PileupSummary(contig, position, refCount, altCount, otherAltCount, alleleFrequency));

        final File file = File.createTempFile("pileup_sumary", ".table");
        PileupSummary.writePileupSummaries(ps, file);
        final List<PileupSummary> psCopy = PileupSummary.readPileupSummaries(file);

        Assert.assertEquals(psCopy.size(), 1);
        Assert.assertEquals(psCopy.get(0).getContig(), contig);
        Assert.assertEquals(psCopy.get(0).getStart(), position);
        Assert.assertEquals(psCopy.get(0).getAltCount(), altCount);
        Assert.assertEquals(psCopy.get(0).getRefCount(), refCount);
        Assert.assertEquals(psCopy.get(0).getOtherAltCount(), otherAltCount);
        Assert.assertEquals(psCopy.get(0).getAlleleFrequency(), alleleFrequency);
    }

}