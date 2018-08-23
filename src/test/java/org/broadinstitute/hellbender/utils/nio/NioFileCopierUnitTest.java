package org.broadinstitute.hellbender.utils.nio;

import org.apache.commons.codec.digest.DigestUtils;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Unit test suite for the {@link NioFileCopier}.
 * Created by jonn on 8/27/18.
 */
public class NioFileCopierUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    private NioFileCopier getNioFileCopierHelper(final Path source, final Path dest, final Boolean overwriteExisting, final Boolean doLogProgress) {
        final NioFileCopier xerox;
        if ( overwriteExisting == null && doLogProgress == null ){
            xerox = NioFileCopier.create(source, dest);
        }
        else if ( doLogProgress == null ) {
            xerox = NioFileCopier.create(source, dest, overwriteExisting);
        }
        else if ( overwriteExisting == null ) {
            xerox = NioFileCopier.create(source, dest, false, doLogProgress);
        }
        else {
            xerox = NioFileCopier.create(source, dest, overwriteExisting, doLogProgress);
        }
        return xerox;
    }

    private Path getSourcePathFromPseudoUrl(final String source) {
        final Path sourcePath;
        if (source.startsWith("GS:")) {
            sourcePath = IOUtils.getPath(getGCPTestInputPath() + source.substring(3));
        }
        else {
            final String sourceBaseName = FilenameUtils.getBaseName(source);
            final String sourceExtension = FilenameUtils.getExtension(source);
            sourcePath = createTempFile(sourceBaseName, sourceExtension).toPath();
        }
        return sourcePath;
    }

    private Path getDestPathFromPseudoUrl(final String dest) {
        final Path destPath;
        if (dest.startsWith("GS:")) {
            final String destBaseName = FilenameUtils.getBaseName(dest.substring(3));
            final String destExtension = FilenameUtils.getExtension(dest.substring(3));
            destPath = BucketUtils.getPathOnGcs(BucketUtils.getTempFilePath(getGCPTestStaging() + destBaseName, destExtension));
        }
        else {
            destPath = getSafeNonExistentPath(dest);
        }
        return destPath;
    }

    private NioFileCopier.ChecksumCalculator createChecksumCalculator(final String checksumCalculatorString) {
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
                    true
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME),
                    createTempFile("tmp2", "fasta").toPath(),
                    false,
                    true
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH),
                    createTempFile("tmp3", "vcf").toPath(),
                    true,
                    false
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.HG38_3_REFERENCE_FILE_NAME),
                    createTempFile("tmp4", "fasta").toPath(),
                    false,
                    false
                },
                {
                    IOUtils.getPath(FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME),
                    createTempFile("tmp5", "vcf").toPath(),
                    null,
                    false
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
                    true,
                    false
            },
            {
                    IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                    createTempFile("tmp", ".tar.gz").toPath(),
                    true,
                    true,
                    true
            },
        };
    }

    @DataProvider
    private Object[][] provideForTestSetDoLogProgress() {

        return new Object[][] {
                {
                        IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                        createTempFile("tmp", ".tar.gz").toPath(),
                        true,
                        false,
                        false
                },
                {
                        IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                        createTempFile("tmp", ".tar.gz").toPath(),
                        true,
                        true,
                        true
                },
                {
                        IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME),
                        createTempFile("tmp", ".tar.gz").toPath(),
                        true,
                        null,
                        true
                },
        };
    }

    @DataProvider
    private Object[][] provideBooleanValuesToSetAndCheckAfterInitialization() {

        final NioFileCopier xerox = getNioFileCopierHelper(IOUtils.getPath(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME), createTempFile("tmp", ".tar.gz").toPath(), true, true);

        return new Object[][]{
                {
                        xerox,
                        false,
                        false
                },
                {
                        xerox,
                        true,
                        true
                },
                {
                        xerox,
                        null,
                        true
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
    private Object[][] provideForTestValidateIntegrity() {

        return new Object[][]{
                // Test multiple checksums:
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        "md5",
                        "325df528d9103d075413e2aa168a8d8a",
                        true
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        "sha256",
                        "f7deb2461c9a576c432262012c06244bb4384ea0fb3eccfc24991932fb1294f2",
                        true
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        "sha512",
                        "406129e1186790b08cb0709d9f6f75f22e3bf13550ab7bc2f11c5f21405edae05955ad83b5de33d28333984a898ddc75b6caad71bd0f58f4a0ca03591ff41479",
                        true
                },
                // Test checksum mismatches / unexpected values:
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "GS:tmp2.dict",
                        "md5",
                        "406129e1186790b08cb0709d9f6f75f22e3bf13550ab7bc2f11c5f21405edae05955ad83b5de33d28333984a898ddc75b6caad71bd0f58f4a0ca03591ff41479",
                        false
                },
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "GS:tmp2.dict",
                        "sha256",
                        "ASDHFJAIWJOF",
                        false
                },
        };
    }

    @DataProvider
    private Object[][] provideForTestGetLatestChecksum() {

        return new Object[][]{
                // Only One Checksum:
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        "md5",
                        "325df528d9103d075413e2aa168a8d8a",
                        null,
                        null
                },
                // Two of the same checksum:
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        "sha256",
                        "f7deb2461c9a576c432262012c06244bb4384ea0fb3eccfc24991932fb1294f2",
                        "sha256",
                        "f7deb2461c9a576c432262012c06244bb4384ea0fb3eccfc24991932fb1294f2",
                },
                // Two different checksums:
                {
                        "GS:benchmark/Solexa-272222.bam.md5",
                        "tmp2.dict",
                        "sha512",
                        "406129e1186790b08cb0709d9f6f75f22e3bf13550ab7bc2f11c5f21405edae05955ad83b5de33d28333984a898ddc75b6caad71bd0f58f4a0ca03591ff41479",
                        "sha256",
                        "f7deb2461c9a576c432262012c06244bb4384ea0fb3eccfc24991932fb1294f2",
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testCreate(final Path source, final Path dest, final Boolean overwriteExisting, final Boolean doLogProgress) {
        final NioFileCopier xerox = getNioFileCopierHelper(source, dest, overwriteExisting, doLogProgress);

        Assert.assertEquals( xerox.getSource(), source );
        Assert.assertEquals( xerox.getDest(), dest );

        if ( overwriteExisting != null ) {
            Assert.assertEquals( xerox.isOverwriteExisting(), overwriteExisting.booleanValue() );
        }
        if ( doLogProgress != null ) {
            Assert.assertEquals( xerox.isDoLogProgress(), doLogProgress.booleanValue() );
        }
    }

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testGetDest(final Path source, final Path dest, final Boolean overwriteExisting, final Boolean doLogProgress) {
        final NioFileCopier xerox = getNioFileCopierHelper(source, dest, overwriteExisting, doLogProgress);

        Assert.assertEquals( xerox.getDest(), dest );
    }

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testGetSource(final Path source, final Path dest, final Boolean overwriteExisting, final Boolean doLogProgress) {
        final NioFileCopier xerox = getNioFileCopierHelper(source, dest, overwriteExisting, doLogProgress);

        Assert.assertEquals( xerox.getSource(), source );
    }

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testIsOverwriteExisting(final Path source, final Path dest, final Boolean overwriteExisting, final Boolean doLogProgress) {
        final NioFileCopier xerox = getNioFileCopierHelper(source, dest, overwriteExisting, doLogProgress);

        if ( overwriteExisting != null ) {
            Assert.assertEquals(xerox.isOverwriteExisting(), overwriteExisting.booleanValue());
        }
        else {
            Assert.assertEquals(xerox.isOverwriteExisting(), false);
        }
    }

    @Test(dataProvider = "provideForTestSetOverwriteExisting")
    public void testSetOverwriteExisting(final Path source, final Path dest, final Boolean overwriteExisting, final Boolean doLogProgress, final boolean expected) {
        final NioFileCopier xerox = getNioFileCopierHelper(source, dest, overwriteExisting, doLogProgress);

        Assert.assertEquals( xerox.isOverwriteExisting(), expected );
    }

    @Test(dataProvider = "provideForTestCreateAndGettersAndSetters")
    public void testIsDoLogProgress(final Path source, final Path dest, final Boolean overwriteExisting, final Boolean doLogProgress) {
        final NioFileCopier xerox = getNioFileCopierHelper(source, dest, overwriteExisting, doLogProgress);

        if ( doLogProgress != null ) {
            Assert.assertEquals(xerox.isDoLogProgress(), doLogProgress.booleanValue());
        }
        else {
            Assert.assertEquals(xerox.isDoLogProgress(), true);
        }
    }

    @Test(dataProvider = "provideForTestSetDoLogProgress")
    public void testSetDoLogProgress(final Path source, final Path dest, final Boolean overwriteExisting, final Boolean doLogProgress, final boolean expected) {
        final NioFileCopier xerox = getNioFileCopierHelper(source, dest, overwriteExisting, doLogProgress);

        Assert.assertEquals( xerox.isDoLogProgress(), expected );
    }

    @Test(dataProvider = "provideBooleanValuesToSetAndCheckAfterInitialization")
    public void testIsLogVerboseGetAndSet(final NioFileCopier xerox, final Boolean testVal, final boolean expected) {
        if ( testVal != null ) {
            xerox.setLogVerbose( testVal );
        }

        Assert.assertEquals( xerox.isLogVerbose(), expected );
    }

    @Test(dataProvider = "provideBooleanValuesToSetAndCheckAfterInitialization")
    public void testIsSilentCopyGetAndSet(final NioFileCopier xerox, final Boolean testVal, final boolean expected) {
        if ( testVal != null ) {
            xerox.setSilentCopy( testVal );
        }

        Assert.assertEquals( xerox.isSilentCopy(), expected );
    }

    @Test(
            dataProvider = "provideForTestInititateCopy",
            groups={"bucket"}
    )
    public void testInitiateCopy(final String source, final String dest) {

        // Get our source:
        final Path sourcePath = getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = getDestPathFromPseudoUrl(dest);

        // Create the copy machine:
        final NioFileCopier xerox = NioFileCopier.create(sourcePath, destPath).setSilentCopy(true);

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
        final Path sourcePath = getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = getDestPathFromPseudoUrl(dest);

        // Copy the file a first time:
        NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).initiateCopy();

        // Do the copy a second time to ensure that dest already exists
        // and set the overwrite flag to ensure it will not fail:
        NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).setOverwriteExisting(true).initiateCopy();

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
        final Path sourcePath = getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = getDestPathFromPseudoUrl(dest);

        // Copy the file a first time:
        NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).initiateCopy();

        // Do the copy a second time to ensure that dest already exists
        // and set the overwrite flag to ensure it will fail:
        NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).setOverwriteExisting(false).initiateCopy();
    }

    @Test(
            dataProvider = "provideForTestInititateCopyOverwrite",
            groups={"bucket"}
    )
    public void testIsCopyComplete(final String source, final String dest) {
        // Get our source:
        final Path sourcePath = getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = getDestPathFromPseudoUrl(dest);

        // Create the xerox machine:
        final NioFileCopier xerox = NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).setOverwriteExisting(true);

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
        final Path sourcePath = getSourcePathFromPseudoUrl("GS:benchmark/Solexa-272222.bam.md5");

        // Get our dest:
        final Path destPath = getDestPathFromPseudoUrl("tmp2.dict");

        // Create the xerox machine:
        final NioFileCopier xerox = NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).setOverwriteExisting(true);

        // Copy the file once:
        xerox.initiateCopy();

        // Attempt to copy the file a second time:
        xerox.initiateCopy();
    }

    @Test(
            dataProvider = "provideForTestValidateIntegrity",
            groups={"bucket"}
    )
    public void testValidateIntegrity(final String source,
                                      final String dest,
                                      final String checksumCalculatorString,
                                      final String expectedChecksum,
                                      final boolean expectedValidity) {
        // Get our source:
        final Path sourcePath = getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = getDestPathFromPseudoUrl(dest);

        // Get our ChecksumCalculator:
        final NioFileCopier.ChecksumCalculator checksumCalculator = createChecksumCalculator(checksumCalculatorString);

        // Create the xerox machine:
        final NioFileCopier xerox = NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).setOverwriteExisting(true);

        // Copy the file:
        xerox.initiateCopy();

        // Validate the integrity of the file:
        Assert.assertEquals(xerox.validateIntegrity(expectedChecksum, checksumCalculator), expectedValidity);
    }

    @Test(
            dataProvider = "provideForTestGetLatestChecksum",
            groups={"bucket"}
    )
    public void testGetLatestChecksum(final String source,
                                      final String dest,
                                      final String checksumCalculatorString,
                                      final String expectedChecksum,
                                      final String secondChecksumCalculatorString,
                                      final String secondExpectedChecksum) {
        // Get our source:
        final Path sourcePath = getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = getDestPathFromPseudoUrl(dest);

        // Get our ChecksumCalculator:
        final NioFileCopier.ChecksumCalculator checksumCalculator = createChecksumCalculator(checksumCalculatorString);

        // Create the xerox machine:
        final NioFileCopier xerox = NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).setOverwriteExisting(true);

        // Copy the file:
        xerox.initiateCopy();

        // Validate the integrity of the file:
        Assert.assertEquals(xerox.validateIntegrity(expectedChecksum, checksumCalculator), true);

        final String lastChecksum;
        if ( secondChecksumCalculatorString != null ) {
            // Validate integrity a second time with the other checksum:

            final NioFileCopier.ChecksumCalculator checksumCalculator2 = createChecksumCalculator(secondChecksumCalculatorString);
            Assert.assertEquals(xerox.validateIntegrity(secondExpectedChecksum, checksumCalculator2), true);

            lastChecksum = secondExpectedChecksum;
        }
        else {
            lastChecksum = expectedChecksum;
        }

        Assert.assertEquals(xerox.getLatestChecksum(), lastChecksum);
    }

    @Test(
            dataProvider = "provideForTestGetLatestChecksum",
            groups={"bucket"}
    )
    public void testGetLatestChecksumCalculator(final String source,
                                                final String dest,
                                                final String checksumCalculatorString,
                                                final String expectedChecksum,
                                                final String secondChecksumCalculatorString,
                                                final String secondExpectedChecksum) {
        // Get our source:
        final Path sourcePath = getSourcePathFromPseudoUrl(source);

        // Get our dest:
        final Path destPath = getDestPathFromPseudoUrl(dest);

        // Get our ChecksumCalculator:
        final NioFileCopier.ChecksumCalculator checksumCalculator = createChecksumCalculator(checksumCalculatorString);

        // Create the xerox machine:
        final NioFileCopier xerox = NioFileCopier.create(sourcePath, destPath).setSilentCopy(true).setOverwriteExisting(true);

        // Copy the file:
        xerox.initiateCopy();

        // Validate the integrity of the file:
        Assert.assertEquals(xerox.validateIntegrity(expectedChecksum, checksumCalculator), true);

        final NioFileCopier.ChecksumCalculator lastChecksumCalculator;
        if ( secondChecksumCalculatorString != null ) {
            // Validate integrity a second time with the other checksum:

            final NioFileCopier.ChecksumCalculator checksumCalculator2 = createChecksumCalculator(secondChecksumCalculatorString);
            Assert.assertEquals(xerox.validateIntegrity(secondExpectedChecksum, checksumCalculator2), true);

            lastChecksumCalculator = checksumCalculator2;
        }
        else {
            lastChecksumCalculator = checksumCalculator;
        }

        Assert.assertEquals(xerox.getLatestChecksumCalculator(), lastChecksumCalculator);
    }
}

