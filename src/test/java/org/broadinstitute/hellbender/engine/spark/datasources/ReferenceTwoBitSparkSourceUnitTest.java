package org.broadinstitute.hellbender.engine.spark.datasources;

import org.broadinstitute.hellbender.engine.GATKPath;
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
        ReferenceSparkSource fastaRef = new ReferenceFileSparkSource(new GATKPath(fastaRefURL));
        ReferenceSparkSource twoBitRef = new ReferenceTwoBitSparkSource(new GATKPath(twoBitRefURL));
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
        final ReferenceTwoBitSparkSource twoBitRef = new ReferenceTwoBitSparkSource(new GATKPath(publicTestDir + "large/human_g1k_v37.20.21.2bit"));
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

    @DataProvider(name="referenceTestCases")
    public Object[][] getReferenceTestCases() {
        return new Object[][] {
                { twoBitRefURL, true },
                { "file://" + twoBitRefURL, true },
                { hg38Reference, false }, // gzipped
                { "file://" + hg38Reference, false }, // gzipped
                { GCS_b37_CHR20_21_REFERENCE_2BIT, true },
                { GCS_b37_CHR20_21_REFERENCE, false },
                // dummy query params at the end to make sure URI.getPath does the right thing
                { GCS_b37_CHR20_21_REFERENCE_2BIT + "?query=param", true },
                { GCS_b37_CHR20_21_REFERENCE + "?query=param", false},
        };
    }

    @Test(dataProvider = "referenceTestCases")
    public void testIsTwoBit(final String referenceSpec, final boolean expectedIsTwoBit) {
        Assert.assertEquals(ReferenceTwoBitSparkSource.isTwoBit(new GATKPath(referenceSpec)), expectedIsTwoBit);
    }

}
