package org.broadinstitute.hellbender.engine.spark.datasources;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;

public class ReferenceTwoBitSparkSourceUnitTest extends GATKBaseTest {
    private static String fastaRefURL = publicTestDir + "large/human_g1k_v37.20.21.fasta";
    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @DataProvider(name = "goodIntervals")
    public Object[][] goodIntervals() throws IOException {
        ReferenceSparkSource fastaRef = new ReferenceFileSparkSource(fastaRefURL);
        ReferenceSparkSource twoBitRef = new ReferenceTwoBitSparkSource(twoBitRefURL);
        return new Object[][]{
                {fastaRef, twoBitRef, "20:2-10"},
                {fastaRef, twoBitRef, "20:4-5"},
                {fastaRef, twoBitRef, "20:4,000-5,000"},
                {fastaRef, twoBitRef, "20:4,0,0,0-5,0,0,0"}, //this is OK too, we just remove commas wherever they are
        };
    }

    @Test(dataProvider = "goodIntervals")
    public void testIntervalConversion( ReferenceSparkSource fastaRef, ReferenceSparkSource twoBitRef, String intervalString) throws IOException {
        SimpleInterval interval = new SimpleInterval(intervalString);
        ReferenceBases expected = fastaRef.getReferenceBases(interval);
        ReferenceBases actual = twoBitRef.getReferenceBases(interval);
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name = "outOfBoundsIntervals")
    public Object[][] getOutOfBoundsIntervals() throws IOException {
        final ReferenceTwoBitSparkSource twoBitRef = new ReferenceTwoBitSparkSource(publicTestDir + "large/human_g1k_v37.20.21.2bit");
        final int chr20End = 63025520;

        return new Object[][] {
                { twoBitRef, new SimpleInterval("20", chr20End - 100, chr20End + 1), 101, chr20End },
                { twoBitRef, new SimpleInterval("20", chr20End - 100, chr20End + 100), 101, chr20End },
                { twoBitRef, new SimpleInterval("20", chr20End - 1, chr20End + 1), 2, chr20End },
                { twoBitRef, new SimpleInterval("20", chr20End, chr20End + 1), 1, chr20End },
        };
    }

    @Test(dataProvider = "outOfBoundsIntervals")
    public void testQueryPastContigEnd( final ReferenceTwoBitSparkSource refSource, final SimpleInterval outOfBoundsInterval, final int expectedNumBases, final int contigEnd ) throws IOException {
        final ReferenceBases bases = refSource.getReferenceBases(outOfBoundsInterval);

        // Verify that the ReferenceTwoBitSource cropped our out-of-bounds interval at the contig end, as expected,
        // and that we got the correct number of bases back.
        Assert.assertEquals(bases.getInterval().getEnd(), contigEnd, "Interval was not cropped at contig end");
        Assert.assertEquals(bases.getBases().length, expectedNumBases, "Wrong number of bases returned from query");
        Assert.assertEquals(bases.getInterval().size(), expectedNumBases, "Wrong interval in ReferenceBases object returned from query");
    }
}
