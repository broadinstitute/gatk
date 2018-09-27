package org.broadinstitute.hellbender.engine.spark.datasources;

import org.apache.spark.SparkConf;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ReferenceMultiSparkSourceUnitTest extends GATKBaseTest {

    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @Test
    public void testSerializeRoundTrip2Bit() {
        ReferenceMultiSparkSource referenceMultiSource = new ReferenceMultiSparkSource(twoBitRefURL, ReferenceWindowFunctions.IDENTITY_FUNCTION);

        final ReferenceMultiSparkSource roundTrippedReference = SparkTestUtils.roundTripInKryo(referenceMultiSource, ReferenceMultiSparkSource.class, new SparkConf());

        Assert.assertEquals(roundTrippedReference.getReferenceSequenceDictionary(null), referenceMultiSource.getReferenceSequenceDictionary(null),
                "\nActual ref: " + roundTrippedReference.getReferenceSequenceDictionary(null) + "\nExpected ref: " + referenceMultiSource.getReferenceSequenceDictionary(null));
        Assert.assertNotNull(roundTrippedReference.getReferenceWindowFunction());
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testBadReferenceFile() {
        new ReferenceMultiSparkSource(
                GATKBaseTest.getSafeNonExistentFile("NonExistentReference.fasta").getAbsolutePath(),
                ReferenceWindowFunctions.IDENTITY_FUNCTION);
    }

}
