package org.broadinstitute.hellbender.engine;

import org.apache.spark.SparkConf;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ReferenceMultiSourceUnitTest extends GATKBaseTest {

    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @Test
    public void testSerializeRoundTrip2Bit() {
        ReferenceMultiSource referenceMultiSource = new ReferenceMultiSource(twoBitRefURL, ReferenceWindowFunctions.IDENTITY_FUNCTION);

        final ReferenceMultiSource roundTrippedReference = SparkTestUtils.roundTripInKryo(referenceMultiSource, ReferenceMultiSource.class, new SparkConf());

        Assert.assertEquals(roundTrippedReference.getReferenceSequenceDictionary(null), referenceMultiSource.getReferenceSequenceDictionary(null),
                "\nActual ref: " + roundTrippedReference.getReferenceSequenceDictionary(null) + "\nExpected ref: " + referenceMultiSource.getReferenceSequenceDictionary(null));
        Assert.assertNotNull(roundTrippedReference.getReferenceWindowFunction());
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testBadReferenceFile() {
        new ReferenceMultiSource(
                GATKBaseTest.getSafeNonExistentFile("NonExistentReference.fasta").getAbsolutePath(),
                ReferenceWindowFunctions.IDENTITY_FUNCTION);
    }

}
