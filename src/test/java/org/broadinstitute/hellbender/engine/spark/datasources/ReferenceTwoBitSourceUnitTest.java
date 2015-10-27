package org.broadinstitute.hellbender.engine.spark.datasources;

import org.broadinstitute.hellbender.engine.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;

public class ReferenceTwoBitSourceUnitTest extends BaseTest {
    private static String fastaRefURL = publicTestDir + "large/human_g1k_v37.20.21.fasta";
    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @DataProvider(name = "goodIntervals")
    public Object[][] goodIntervals() throws IOException {
        ReferenceSource fastaRef = new ReferenceFileSource(fastaRefURL);
        ReferenceSource twoBitRef = new ReferenceTwoBitSource(null, twoBitRefURL);
        return new Object[][]{
                {fastaRef, twoBitRef, "20:2-10"},
                {fastaRef, twoBitRef, "20:4-5"},
                {fastaRef, twoBitRef, "20:4,000-5,000"},
                {fastaRef, twoBitRef, "20:4,0,0,0-5,0,0,0"}, //this is OK too, we just remove commas wherever they are
        };
    }

    @Test(dataProvider = "goodIntervals")
    public void testIntervalConversion(ReferenceSource fastaRef, ReferenceSource twoBitRef, String intervalString) throws IOException {
        SimpleInterval interval = new SimpleInterval(intervalString);
        ReferenceBases expected = fastaRef.getReferenceBases(null, interval);
        ReferenceBases actual = twoBitRef.getReferenceBases(null, interval);
        Assert.assertEquals(actual, expected);
    }
}
