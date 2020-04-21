package org.broadinstitute.hellbender.engine.spark.datasources;

import org.apache.spark.SparkConf;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class ReferenceMultiSparkSourceUnitTest extends GATKBaseTest {
    private String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @DataProvider(name="referenceTestCases")
    public Object[][] getReferenceTestCases() {
        return new Object[][] {
                { twoBitRefURL, false },
                { "file://" + twoBitRefURL, false },
                { hg38Reference, true }, // gzipped
                { "file://" + hg38Reference, true }, // gzipped
                { GCS_b37_CHR20_21_REFERENCE_2BIT, false },
                { GCS_b37_CHR20_21_REFERENCE, true },
                // dummy query params at the end to make sure URI.getPath does the right thing
                { GCS_b37_CHR20_21_REFERENCE + "?query=param", true}
        };
    }

    @Test(dataProvider = "referenceTestCases")
    public void testIsFasta(final String referenceSpec, final boolean expectedIsFasta) {
        Assert.assertEquals(ReferenceMultiSparkSource.isFasta(new GATKPathSpecifier(referenceSpec)), expectedIsFasta);
    }

    @Test
    public void testSerializeRoundTrip2Bit() {
        ReferenceMultiSparkSource referenceMultiSource = new ReferenceMultiSparkSource(new GATKPathSpecifier(twoBitRefURL), ReferenceWindowFunctions.IDENTITY_FUNCTION);

        final ReferenceMultiSparkSource roundTrippedReference = SparkTestUtils.roundTripInKryo(referenceMultiSource, ReferenceMultiSparkSource.class, new SparkConf());

        Assert.assertEquals(roundTrippedReference.getReferenceSequenceDictionary(null), referenceMultiSource.getReferenceSequenceDictionary(null),
                "\nActual ref: " + roundTrippedReference.getReferenceSequenceDictionary(null) + "\nExpected ref: " + referenceMultiSource.getReferenceSequenceDictionary(null));
        Assert.assertNotNull(roundTrippedReference.getReferenceWindowFunction());
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testBadReferenceFile() {
        new ReferenceMultiSparkSource(
                new GATKPathSpecifier(GATKBaseTest.getSafeNonExistentFile("NonExistentReference.fasta").getAbsolutePath()),
                ReferenceWindowFunctions.IDENTITY_FUNCTION);
    }

}
