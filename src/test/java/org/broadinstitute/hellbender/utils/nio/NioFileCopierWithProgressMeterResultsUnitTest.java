package org.broadinstitute.hellbender.utils.nio;

import org.apache.commons.codec.digest.MessageDigestAlgorithms;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;

/**
 * Unit test suite for {@link NioFileCopierWithProgressMeterResults}.
 * Created by jonn on 9/4/18.
 */
public class NioFileCopierWithProgressMeterResultsUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Data Providers:

    @DataProvider private Object[][] provideForTestFieldGetters() {
        return new Object[][] {
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        32,
                        MessageDigestAlgorithms.MD5,
                        "325df528d9103d075413e2aa168a8d8a",
                        true
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        32,
                        MessageDigestAlgorithms.SHA_256,
                        "f7deb2461c9a576c432262012c06244bb4384ea0fb3eccfc24991932fb1294f2",
                        true
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "GS:tmp2.dict",
                        32,
                        MessageDigestAlgorithms.MD5,
                        "406129e1186790b08cb0709d9f6f75f22e3bf13550ab7bc2f11c5f21405edae05955ad83b5de33d28333984a898ddc75b6caad71bd0f58f4a0ca03591ff41479",
                        false
                },
                {
                        "GS:benchmark/human_g1k_v37.dict",
                        "GS:testOut.dict",
                        10005,
                        MessageDigestAlgorithms.MD5,
                        "8f04f9288bd8876b8b76a83917ac2725",
                        true
                }
        };
    }

    @DataProvider private Object[][] provideForTestFieldGettersNoValidation() {
        return new Object[][] {
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        32,
                },
                {
                        "GS:benchmark/human_g1k_v37.dict",
                        "GS:testOut.dict",
                        10005,
                }
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(
            dataProvider = "provideForTestFieldGetters",
            groups={"bucket"}
    )
    public void testFieldGetters(final String source,
                                  final String dest,
                                  final long expectedSize,
                                  final String checksumAlgorithm,
                                  final String expectedChecksum,
                                  final boolean expectedValidity) {
        // Get our source:
        final Path sourcePath = NioFileCopierWithProgessMeterTestUtils.getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = NioFileCopierWithProgessMeterTestUtils.getDestPathFromPseudoUrl(dest);

        // Create the xerox machine:
        final NioFileCopierWithProgressMeter xerox =
                NioFileCopierWithProgressMeter.create(sourcePath, destPath)
                        .setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT)
                        .setOverwriteExisting(true)
                        .setChecksumAlgorithmAndExpectedChecksum(checksumAlgorithm, expectedChecksum);

        // Copy the file:
        final NioFileCopierWithProgressMeterResults results = xerox.initiateCopy();

        // Validate that all the getters have the correct values:
        Assert.assertEquals( results.getSource(), sourcePath );
        Assert.assertEquals( results.getDestination(), destPath );
        Assert.assertEquals( results.getSize(), expectedSize );
        Assert.assertEquals( results.wasValidationRequested(), true );
        Assert.assertEquals( results.getChecksumAlgorithm(), checksumAlgorithm );
        Assert.assertEquals( results.getExpectedChecksum(), expectedChecksum );
        Assert.assertEquals( results.isDestFileValid(), expectedValidity );
    }

    @Test(
            dataProvider = "provideForTestFieldGettersNoValidation",
            groups={"bucket"}
    )
    public void testFieldGettersNoValidation(final String source,
                                 final String dest,
                                 final long expectedSize) {
        // Get our source:
        final Path sourcePath = NioFileCopierWithProgessMeterTestUtils.getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = NioFileCopierWithProgessMeterTestUtils.getDestPathFromPseudoUrl(dest);

        // Create the xerox machine:
        final NioFileCopierWithProgressMeter xerox =
                NioFileCopierWithProgressMeter.create(sourcePath, destPath)
                        .setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT)
                        .setOverwriteExisting(true);

        // Copy the file:
        final NioFileCopierWithProgressMeterResults results = xerox.initiateCopy();

        // Validate that all the getters have the correct values:
        Assert.assertEquals( results.getSource(), sourcePath );
        Assert.assertEquals( results.getDestination(), destPath );
        Assert.assertEquals( results.getSize(), expectedSize );
        Assert.assertEquals( results.wasValidationRequested(), false );
        Assert.assertEquals( results.getChecksumAlgorithm(), "" );
        Assert.assertEquals( results.getExpectedChecksum(), "" );
        Assert.assertEquals( results.isDestFileValid(), false );
    }

}
