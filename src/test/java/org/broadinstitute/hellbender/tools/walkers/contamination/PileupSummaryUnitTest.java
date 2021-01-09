package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Created by David Benjamin on 2/15/17.
 */
public class PileupSummaryUnitTest extends GATKBaseTest {

    @Test
    public void test() {
        final String contig = "chr1";
        final int position = 1;
        final int refCount = 20;
        final int altCount = 10;
        final int otherAltCount = 2;
        final double alleleFrequency = 0.3;

        final List<PileupSummary> ps = Arrays.asList(new PileupSummary(contig, position, refCount, altCount, otherAltCount, alleleFrequency));

        final File file = IOUtils.createTempFile("pileup_sumary", ".table");
        PileupSummary.writeToFile("sample", ps, file);
        final List<PileupSummary> psCopy = PileupSummary.readFromFile(file).getRight();

        Assert.assertEquals(psCopy.size(), 1);
        Assert.assertEquals(psCopy.get(0).getContig(), contig);
        Assert.assertEquals(psCopy.get(0).getStart(), position);
        Assert.assertEquals(psCopy.get(0).getAltCount(), altCount);
        Assert.assertEquals(psCopy.get(0).getRefCount(), refCount);
        Assert.assertEquals(psCopy.get(0).getOtherAltCount(), otherAltCount);
        Assert.assertEquals(psCopy.get(0).getAlleleFrequency(), alleleFrequency);
    }

    @DataProvider
    public Object[][] comparatorData(){
        return new Object[][]{
                // contig1, position1, contig2, position2, expected
                { "MT", 3,   "MT", 3,   0 },
                { "1",  3,   "2",  2,   -1 },
                { "22", 3,   "X",  2,   -1 },
                { "21", 3,   "MT", 2,   -1 },
                { "MT", 3,   "X",  2,   1 },
                { "X",  300, "X",  800, -1 },
        };
    }

    @Test(dataProvider = "comparatorData")
    public void testComparator(final String contig1, final int position1,
                               final String contig2, final int position2,
                               int expected){
        final SAMSequenceDictionary dict = new ReferenceFileSource(new File(b37Reference).toPath()).getSequenceDictionary();
        final PileupSummary ps1 = new PileupSummary(contig1, position1, 0,0,0, 0);
        final PileupSummary ps2 = new PileupSummary(contig2, position2, 0,0,0, 0);

        // They must both be positive or negative - so use exlusive-OR
        Assert.assertFalse(new PileupSummary.PileupSummaryComparator(dict).compare(ps1, ps2) < 0 ^ expected < 0);
    }
}