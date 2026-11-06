package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
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
        Assert.assertEquals(psCopy.get(0).getRefBase(), Nucleotide.N);
        Assert.assertEquals(psCopy.get(0).getAltBase(), Nucleotide.N);
    }

    @Test
    public void testBases() {
        final List<PileupSummary> ps = Arrays.asList(new PileupSummary("chr1", 1, 20, 10, 2, 0.3, Nucleotide.C, Nucleotide.T));

        final File file = IOUtils.createTempFile("pileup_summary", ".table");
        PileupSummary.writeToFile("sample", ps, file);
        final PileupSummary psCopy = PileupSummary.readFromFile(file).getRight().get(0);

        Assert.assertEquals(psCopy.getRefBase(), Nucleotide.C);
        Assert.assertEquals(psCopy.getAltBase(), Nucleotide.T);
    }

    // tables written before the bases were recorded can still be read
    @Test
    public void testTableWithoutBases() throws IOException {
        final File file = IOUtils.createTempFile("pileup_summary", ".table");
        Files.write(file.toPath(), Arrays.asList("#<METADATA>SAMPLE=sample",
                "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency",
                "chr1\t1\t20\t10\t2\t0.3"));
        final PileupSummary ps = PileupSummary.readFromFile(file).getRight().get(0);

        Assert.assertEquals(ps.getStart(), 1);
        Assert.assertEquals(ps.getRefCount(), 20);
        Assert.assertEquals(ps.getAltCount(), 10);
        Assert.assertEquals(ps.getOtherAltCount(), 2);
        Assert.assertEquals(ps.getAlleleFrequency(), 0.3);
        Assert.assertEquals(ps.getRefBase(), Nucleotide.N);
        Assert.assertEquals(ps.getAltBase(), Nucleotide.N);
    }

    @Test
    public void testReadWithFilter() {
        final List<PileupSummary> ps = Arrays.asList(
                new PileupSummary("chr1", 1, 20, 10, 0, 0.1),
                new PileupSummary("chr1", 2, 20, 10, 0, 0.3),
                new PileupSummary("chr1", 3, 20, 10, 0, 0.05));

        final File file = IOUtils.createTempFile("pileup_summary", ".table");
        PileupSummary.writeToFile("sample", ps, file);
        final List<PileupSummary> filtered = PileupSummary.readFromFile(file, site -> site.getAlleleFrequency() < 0.2).getRight();

        Assert.assertEquals(filtered.stream().mapToInt(PileupSummary::getStart).toArray(), new int[] {1, 3});
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