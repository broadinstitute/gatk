package org.broadinstitute.hellbender.utils.nio;

import org.apache.commons.codec.digest.DigestUtils;
import org.apache.commons.codec.digest.MessageDigestAlgorithms;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Unit test suite for the {@link NioFileCopierWithProgressMeter}.
 * Created by jonn on 8/27/18.
 */
public class NioFileCopierWithProgressMeterUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Helper Methods:

    private NioFileCopierWithProgressMeter getNioFileCopierHelper(final Path source, final Path dest, final Boolean overwriteExisting, final NioFileCopierWithProgressMeter.Verbosity verbosity) {

        final NioFileCopierWithProgressMeter xerox;
        if ( overwriteExisting == null && verbosity == null ){
            xerox = NioFileCopierWithProgressMeter.create(source, dest);
        }
        else if ( verbosity == null ) {
            xerox = NioFileCopierWithProgressMeter.create(source, dest, overwriteExisting);
        }
        else if ( overwriteExisting == null ) {
            xerox = NioFileCopierWithProgressMeter.create(source, dest, false, verbosity);
        }
        else {
            xerox = NioFileCopierWithProgressMeter.create(source, dest, overwriteExisting, verbosity);
        }
        return xerox;
    }

    private NioFileCopierWithProgressMeter.ChecksumCalculator createChecksumCalculator(final String checksumCalculatorString) {
        switch (checksumCalculatorString) {
            case "md2":    return DigestUtils::md2Hex;
            case "md5":    return DigestUtils::md5Hex;
            case "sha1":   return DigestUtils::sha1Hex;
            case "sha256": return DigestUtils::sha256Hex;
            case "sha384": return DigestUtils::sha384Hex;
            case "sha512": return DigestUtils::sha512Hex;
            default:       return DigestUtils::sha256Hex;
        }
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestCreateAndGettersAndSetters() {
        return new Object[][] {
                {
                    IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                    createTempFile("tmp", ".tar.gz").toPath(),
                    true,
                    NioFileCopierWithProgressMeter.Verbosity.MODERATE
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME),
                    createTempFile("tmp2", "fasta").toPath(),
                    false,
                    NioFileCopierWithProgressMeter.Verbosity.MODERATE
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH),
                    createTempFile("tmp3", "vcf").toPath(),
                    true,
                    NioFileCopierWithProgressMeter.Verbosity.MINIMAL
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.HG38_3_REFERENCE_FILE_NAME),
                    createTempFile("tmp4", "fasta").toPath(),
                    false,
                    NioFileCopierWithProgressMeter.Verbosity.MINIMAL
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME),
                    createTempFile("tmp5", "vcf").toPath(),
                    null,
                    NioFileCopierWithProgressMeter.Verbosity.MINIMAL
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.XSV_DEADBEEFSV_FILE_PATH),
                    createTempFile("tmp6", "xsv").toPath(),
                    false,
                    null
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19),
                    createTempFile("tmp7", "vcf").toPath(),
                    null,
                    null
                },
        };
    }

    @DataProvider
    private Object[][] provideForTestSetOverwriteExisting() {

        return new Object[][] {
            {
                    IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                    createTempFile("tmp", ".tar.gz").toPath(),
                    false,
                    NioFileCopierWithProgressMeter.Verbosity.MODERATE,
                    false
            },
            {
                    IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                    createTempFile("tmp", ".tar.gz").toPath(),
                    true,
                    NioFileCopierWithProgressMeter.Verbosity.MODERATE,
                    true
            },
        };
    }

    @DataProvider
    private Object[][] provideForTestSetVerbosity() {

        return new Object[][] {
                {
                        IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                        createTempFile("tmp", ".tar.gz").toPath(),
                        true,
                        NioFileCopierWithProgressMeter.Verbosity.MODERATE,
                        NioFileCopierWithProgressMeter.Verbosity.MODERATE
                },
                {
                        IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                        createTempFile("tmp", ".tar.gz").toPath(),
                        true,
                        NioFileCopierWithProgressMeter.Verbosity.MINIMAL,
                        NioFileCopierWithProgressMeter.Verbosity.MINIMAL
                },
                {
                        IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                        createTempFile("tmp", ".tar.gz").toPath(),
                        true,
                        NioFileCopierWithProgressMeter.Verbosity.SILENT,
                        NioFileCopierWithProgressMeter.Verbosity.SILENT
                },
                {
                        IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                        createTempFile("tmp", ".tar.gz").toPath(),
                        true,
                        NioFileCopierWithProgressMeter.Verbosity.VERBOSE,
                        NioFileCopierWithProgressMeter.Verbosity.VERBOSE
                },
                {
                        IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                        createTempFile("tmp", ".tar.gz").toPath(),
                        true,
                        null,
                        NioFileCopierWithProgressMeter.Verbosity.MODERATE
                },
        };
    }

    @DataProvider
    private Object[][] provideForTestInititateCopy() {

        return new Object[][]{
                {
                        "GS:benchmark/human_g1k_v37.dict",
                        "tmp.dict"
                },
                {
                        "GS:benchmark/human_g1k_v37.dict",
                        "GS:testOut.dict"
                },
                {
                        FuncotatorTestConstants.XSV_DEADBEEFSV_FILE_PATH,
                        "GS:testOut.dict"
                },
                {
                        FuncotatorTestConstants.XSV_DEADBEEFSV_FILE_PATH,
                        "test.xsv"
                },
        };
    }

    @DataProvider
    private Object[][] provideForTestInititateCopyOverwrite() {

        return new Object[][]{
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict"
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "GS:testOut2.dict"
                },
                {
                        FuncotatorTestConstants.XSV_DEADBEEFSV_FILE_PATH,
                        "GS:testOut2.dict"
                },
                {
                        FuncotatorTestConstants.XSV_DEADBEEFSV_FILE_PATH,
                        "test2.xsv"
                },
        };
    }

    @DataProvider
    private Object[][] provideForTestChecksumsAndValidateIntegrity() {

        return new Object[][]{
                // Test multiple checksums:
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        MessageDigestAlgorithms.MD5,
                        "325df528d9103d075413e2aa168a8d8a",
                        true
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        MessageDigestAlgorithms.SHA_256,
                        "f7deb2461c9a576c432262012c06244bb4384ea0fb3eccfc24991932fb1294f2",
                        true
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        MessageDigestAlgorithms.SHA_512,
                        "406129e1186790b08cb0709d9f6f75f22e3bf13550ab7bc2f11c5f21405edae05955ad83b5de33d28333984a898ddc75b6caad71bd0f58f4a0ca03591ff41479",
                        true
                },
                // Test checksum mismatches / unexpected values:
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "GS:tmp2.dict",
                        MessageDigestAlgorithms.MD5,
                        "406129e1186790b08cb0709d9f6f75f22e3bf13550ab7bc2f11c5f21405edae05955ad83b5de33d28333984a898ddc75b6caad71bd0f58f4a0ca03591ff41479",
                        false
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "GS:tmp2.dict",
                        MessageDigestAlgorithms.SHA_256,
                        "ASDHFJAIWJOF",
                        false
                },
        };
    }

    @DataProvider
    private Object[][] provideForTestChecksumsAndValidateIntegrityNoAlgorithm() {

        return new Object[][]{
                // Test multiple checksums:
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testCreate(final Path source, final Path dest, final Boolean overwriteExisting, final NioFileCopierWithProgressMeter.Verbosity verbosity) {
        final NioFileCopierWithProgressMeter xerox = getNioFileCopierHelper(source, dest, overwriteExisting, verbosity);

        Assert.assertEquals( xerox.getSource(), source );
        Assert.assertEquals( xerox.getDest(), dest );

        if ( overwriteExisting != null ) {
            Assert.assertEquals( xerox.isOverwriteExisting(), overwriteExisting.booleanValue() );
        }
        if ( verbosity != null ) {
            Assert.assertEquals( xerox.getVerbosity(), verbosity );
        }
    }

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testGetDest(final Path source, final Path dest, final Boolean overwriteExisting, final NioFileCopierWithProgressMeter.Verbosity verbosity) {
        final NioFileCopierWithProgressMeter xerox = getNioFileCopierHelper(source, dest, overwriteExisting, verbosity);

        Assert.assertEquals( xerox.getDest(), dest );
    }

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testGetSource(final Path source, final Path dest, final Boolean overwriteExisting, final NioFileCopierWithProgressMeter.Verbosity verbosity) {
        final NioFileCopierWithProgressMeter xerox = getNioFileCopierHelper(source, dest, overwriteExisting, verbosity);

        Assert.assertEquals( xerox.getSource(), source );
    }

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testIsOverwriteExisting(final Path source, final Path dest, final Boolean overwriteExisting, final NioFileCopierWithProgressMeter.Verbosity verbosity) {
        final NioFileCopierWithProgressMeter xerox = getNioFileCopierHelper(source, dest, overwriteExisting, verbosity);

        if ( overwriteExisting != null ) {
            Assert.assertEquals(xerox.isOverwriteExisting(), overwriteExisting.booleanValue());
        }
        else {
            Assert.assertEquals(xerox.isOverwriteExisting(), false);
        }
    }

    @Test(dataProvider = "provideForTestSetOverwriteExisting")
    public void testSetOverwriteExisting(final Path source, final Path dest, final Boolean overwriteExisting, final NioFileCopierWithProgressMeter.Verbosity verbosity, final boolean expected) {
        final NioFileCopierWithProgressMeter xerox = getNioFileCopierHelper(source, dest, overwriteExisting, verbosity);

        Assert.assertEquals( xerox.isOverwriteExisting(), expected );
    }

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testIsDoLogProgress(final Path source, final Path dest, final Boolean overwriteExisting, final NioFileCopierWithProgressMeter.Verbosity verbosity) {
        final NioFileCopierWithProgressMeter xerox = getNioFileCopierHelper(source, dest, overwriteExisting, verbosity);

        if ( verbosity != null ) {
            Assert.assertEquals(xerox.getVerbosity(), verbosity);
        }
        else {
            Assert.assertEquals(xerox.getVerbosity(), NioFileCopierWithProgressMeter.Verbosity.MODERATE);
        }
    }

    @Test(dataProvider = "provideForTestSetVerbosity")
    public void testSetVerbosity(final Path source, final Path dest, final Boolean overwriteExisting, final NioFileCopierWithProgressMeter.Verbosity verbosity, final NioFileCopierWithProgressMeter.Verbosity expected) {
        final NioFileCopierWithProgressMeter xerox = getNioFileCopierHelper(source, dest, overwriteExisting, verbosity);

        Assert.assertEquals( xerox.getVerbosity(), expected );
    }

    @Test(
            dataProvider = "provideForTestInititateCopy",
            groups={"bucket"}
    )
    public void testInitiateCopy(final String source, final String dest) {

        // Get our source:
        final Path sourcePath = NioFileCopierWithProgessMeterTestUtils.getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = NioFileCopierWithProgessMeterTestUtils.getDestPathFromPseudoUrl(dest);

        // Create the copy machine:
        final NioFileCopierWithProgressMeter xerox = NioFileCopierWithProgressMeter.create(sourcePath, destPath).setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT);

        // Do the copy:
        xerox.initiateCopy();

        // Verify the output file exists and that the files are the same:
        Assert.assertTrue(Files.exists(destPath));
        try {
            Assert.assertEquals(Files.readAllBytes(destPath), Files.readAllBytes(sourcePath));
        }
        catch (final IOException ex) {
            throw new GATKException("FAILURE: File contents could not be verified: " + destPath.toUri() + " -> " + sourcePath.toUri());
        }
    }

    @Test(
            dataProvider = "provideForTestInititateCopyOverwrite",
            groups={"bucket"}
    )
    public void testInitiateCopyFileExistsOk(final String source, final String dest) {

        // Get our source:
        final Path sourcePath = NioFileCopierWithProgessMeterTestUtils.getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = NioFileCopierWithProgessMeterTestUtils.getDestPathFromPseudoUrl(dest);

        // Copy the file a first time:
        NioFileCopierWithProgressMeter.create(sourcePath, destPath).setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT).initiateCopy();

        // Do the copy a second time to ensure that dest already exists
        // and set the overwrite flag to ensure it will not fail:
        NioFileCopierWithProgressMeter.create(sourcePath, destPath).setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT).setOverwriteExisting(true).initiateCopy();

        // Verify the output file exists and that the files are the same:
        Assert.assertTrue(Files.exists(destPath));
        try {
            Assert.assertEquals(Files.readAllBytes(destPath), Files.readAllBytes(sourcePath));
        }
        catch (final IOException ex) {
            throw new GATKException("FAILURE: File contents could not be verified: " + destPath.toUri() + " -> " + sourcePath.toUri());
        }
    }

    @Test(
            dataProvider = "provideForTestInititateCopyOverwrite",
            groups={"bucket"},
            expectedExceptions = { UserException.CouldNotCreateOutputFile.class}
    )
    public void testInitiateCopyFileExistsException(final String source, final String dest) {

        // Get our source:
        final Path sourcePath = NioFileCopierWithProgessMeterTestUtils.getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = NioFileCopierWithProgessMeterTestUtils.getDestPathFromPseudoUrl(dest);

        // Copy the file a first time:
        NioFileCopierWithProgressMeter.create(sourcePath, destPath).setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT).initiateCopy();

        // Do the copy a second time to ensure that dest already exists
        // and set the overwrite flag to ensure it will fail:
        NioFileCopierWithProgressMeter.create(sourcePath, destPath).setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT).setOverwriteExisting(false).initiateCopy();
    }

    @Test(
            dataProvider = "provideForTestInititateCopyOverwrite",
            groups={"bucket"}
    )
    public void testIsCopyComplete(final String source, final String dest) {

        // Get our source:
        final Path sourcePath = NioFileCopierWithProgessMeterTestUtils.getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = NioFileCopierWithProgessMeterTestUtils.getDestPathFromPseudoUrl(dest);

        // Create the xerox machine:
        final NioFileCopierWithProgressMeter xerox = NioFileCopierWithProgressMeter.create(sourcePath, destPath).setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT).setOverwriteExisting(true);

        Assert.assertFalse( xerox.isCopyComplete() );

        // Copy the file:
        xerox.initiateCopy();

        Assert.assertTrue( xerox.isCopyComplete() );
    }

    @Test(
            groups={"bucket"},
            expectedExceptions = {GATKException.class}
    )
    public void testCanOnlyCopyOnce() {

        // Get our source:
        final Path sourcePath = NioFileCopierWithProgessMeterTestUtils.getSourcePathFromPseudoUrl("GS:benchmark/Solexa-272222.bam.md5");

        // Get our dest:
        final Path destPath = NioFileCopierWithProgessMeterTestUtils.getDestPathFromPseudoUrl("tmp2.dict");

        // Create the xerox machine:
        final NioFileCopierWithProgressMeter xerox = NioFileCopierWithProgressMeter.create(sourcePath, destPath).setVerbosity(NioFileCopierWithProgressMeter.Verbosity.SILENT).setOverwriteExisting(true);

        // Copy the file once:
        xerox.initiateCopy();

        // Attempt to copy the file a second time:
        xerox.initiateCopy();
    }

    @Test(
            dataProvider = "provideForTestChecksumsAndValidateIntegrity",
            groups={"bucket"}
    )
    public void testChecksumsAndValidateIntegrity(final String source,
                                      final String dest,
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

        // Validate the checksum algorithm was set:
        Assert.assertEquals( results.getChecksumAlgorithm(), checksumAlgorithm );

        // Validate the calculated validity of the file:
        Assert.assertEquals( results.isDestFileValid(), expectedValidity );

        // Validate the expected checksum value was set:
        Assert.assertEquals( results.getExpectedChecksum(), expectedChecksum );
    }

    @Test(
            dataProvider = "provideForTestChecksumsAndValidateIntegrityNoAlgorithm",
            groups={"bucket"}
    )
    public void testChecksumsAndValidateIntegrityNoAlgorithm(final String source,
                                                             final String dest) {
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

        // Since we did not set the algorithm to validate the checksum, assert that all the associated values are the
        // "empty" values:
        Assert.assertEquals( results.wasValidationRequested(), false );
        Assert.assertEquals( results.getChecksumAlgorithm(), "" );
        Assert.assertEquals( results.isDestFileValid(), false);
        Assert.assertEquals( results.getExpectedChecksum(), "" );
    }
}

