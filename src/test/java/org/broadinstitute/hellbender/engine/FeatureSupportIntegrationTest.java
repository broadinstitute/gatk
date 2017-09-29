package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithVariants;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystemException;
import java.time.Period;
import java.time.ZonedDateTime;
import java.util.Arrays;

public final class FeatureSupportIntegrationTest extends CommandLineProgramTest {
    private static final String FEATURE_INTEGRATION_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @Override
    public String getTestedClassName() {
        return ExampleReadWalkerWithVariants.class.getSimpleName();
    }

    /**
     * Updates the given file's modification time to the current time.
     * @param f File in which to change the modification time.
     * @return true if the modification time was set; false otherwise.
     */
    public boolean updateFileModifiedTime(final File f) {
        return updateFileModifiedTime(f, ZonedDateTime.now());
    }

    /**
     * Updates the given file's modification time.
     * @param f File in which to change the modification time.
     * @param time Time to which to set the modification time of {@code f}
     * @return true if the modification time was set; false otherwise.
     */
    public boolean updateFileModifiedTime(final File f, final ZonedDateTime time) {
        return f.setLastModified(time.toInstant().toEpochMilli());
    }

    /**
     * Sets the given input file's modification time to be before the given index file's modification time.
     * @param inputFile Input file of which to change the modification time.
     * @param indexFile Index file of which to change the modification time (to some time after that of {@code inputFile}.
     */
    public void testCheckIndexModificationTimeHelper(final File inputFile, final File indexFile) throws FileSystemException {
        if (! updateFileModifiedTime(inputFile) ) {
            throw new FileSystemException("Could not change the time of the given input file: " + inputFile.getAbsolutePath());
        }

        final ZonedDateTime t = ZonedDateTime.now().minus(Period.ofDays(15));
        if (! updateFileModifiedTime(indexFile, t) ) {
            throw new FileSystemException("Could not change the time of the given input file: " + inputFile.getAbsolutePath());
        }
    }

    @Test
    public void testFeatureSupportUsingVCF() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -V " + "TestFeatures:" + FEATURE_INTEGRATION_TEST_DIRECTORY + "feature_data_source_test.vcf" +
                " --groupVariantsBySource" +
                " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeatureSupportUsingVCF_output.txt")
        );
        testSpec.executeTest("testFeatureSupportUsingVCF", this);
    }

    @Test
    public void testFeaturesAsIntervals() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -L " + publicTestDir + "org/broadinstitute/hellbender/utils/interval/intervals_from_features_test.vcf" +
                " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeaturesAsIntervals_output.txt")
        );
        testSpec.executeTest("testFeaturesAsIntervals", this);
    }

    @Test
    public void testFeaturesAsIntervalsWithExclusion() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -L " + publicTestDir + "org/broadinstitute/hellbender/utils/interval/intervals_from_features_test.vcf" +
                " -XL " + publicTestDir + "org/broadinstitute/hellbender/utils/interval/intervals_from_features_test_exclude.vcf" +
                " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeaturesAsIntervalsWithExclusion_output.txt")
        );
        testSpec.executeTest("testFeaturesAsIntervals", this);
    }

    @Test
    public void testFeaturesAsIntervalsNonExistentFile() throws IOException {
        // Non-existent files should be interpreted as interval strings and fail interval parsing
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -L " + publicTestDir + "non_existent_file.vcf" +
                " -O %s",
                1,
                UserException.MalformedGenomeLoc.class
        );
        testSpec.executeTest("testFeaturesAsIntervals", this);
    }

    @Test
    public void testFeaturesAsIntervalsUnrecognizedFormatFile() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                " -L " + publicTestDir + "org/broadinstitute/hellbender/utils/interval/unrecognized_format_file.xyz" +
                " -O %s",
                1,
                UserException.CouldNotReadInputFile.class
        );
        testSpec.executeTest("testFeaturesAsIntervals", this);
    }

    @Test(groups = {"bucket"})
    public void testFeatureInputWithOutOfDateIndexGcsVcf() throws IOException {

        // Non-existent files should be interpreted as interval strings and fail interval parsing
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                        " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                        " -V " + "TestFeatures:" + GCS_GATK_TEST_RESOURCES + "large/emptyDatedFile.vcf" +
                        " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeatureInputWithOutOfDateIndexVcf_output.txt")
        );
        testSpec.executeTest("testVcfIndexOutOfDate", this);
    }

    @Test(groups = {"bucket"})
    public void testFeatureInputWithOutOfDateIndexGcsVcfWithException() throws IOException {

        // Non-existent files should be interpreted as interval strings and fail interval parsing
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                        " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                        " -V " + "TestFeatures:" + GCS_GATK_TEST_RESOURCES + "large/emptyDatedFile.vcf" +
                        " --errorOnOutOfDateIndex true" +
                        " -O %s",
                1,
                UserException.OutOfDateIndex.class
        );
        testSpec.executeTest("testVcfIndexOutOfDate", this);
    }

    @Test()
    public void testFeatureInputWithOutOfDateIndexVcf() throws IOException {

        final File vcf = new File(publicTestDir, "emptyDatedFile.vcf");
        final File vcfIdx = new File(publicTestDir, "emptyDatedFile.vcf.idx");

        // Set the modification time of vcfIdx to be before vcf:
        testCheckIndexModificationTimeHelper(vcf, vcfIdx);

        // Non-existent files should be interpreted as interval strings and fail interval parsing
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                        " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                        " -V " + "TestFeatures:" + vcf.getAbsolutePath() +
                        " -O %s",
                Arrays.asList(FEATURE_INTEGRATION_TEST_DIRECTORY + "expected_testFeatureInputWithOutOfDateIndexVcf_output.txt")
        );
        testSpec.executeTest("testVcfIndexOutOfDate", this);
    }

    @Test()
    public void testFeatureInputWithOutOfDateIndexVcfWithException() throws IOException {

        final File vcf = new File(publicTestDir, "emptyDatedFile.vcf");
        final File vcfIdx = new File(publicTestDir, "emptyDatedFile.vcf.idx");

        // Set the modification time of vcfIdx to be before vcf:
        testCheckIndexModificationTimeHelper(vcf, vcfIdx);

        // Non-existent files should be interpreted as interval strings and fail interval parsing
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg19MiniReference +
                        " -I " + FEATURE_INTEGRATION_TEST_DIRECTORY + "reads_data_source_test1.bam" +
                        " -V " + "TestFeatures:" + vcf.getAbsolutePath() +
                        " --errorOnOutOfDateIndex true" +
                        " -O %s",
                1,
                UserException.OutOfDateIndex.class
        );
        testSpec.executeTest("testVcfIndexOutOfDate", this);
    }
}
