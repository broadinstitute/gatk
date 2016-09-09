package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class TargetCoverageSexGenotyperIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/sexgenotyper/";
    private static final File TEST_CONTIG_PLOIDY_ANNOTS_FILE = new File(TEST_SUB_DIR, "contig_annots.tsv");
    private static final File TEST_RCC_FILE = new File(TEST_SUB_DIR, "sex_genotyper_rcc_trunc.tsv");
    private static final File TEST_SEX_GENOTYPE_FILE = new File(TEST_SUB_DIR, "sex_genotypes_agilent_trunc.tsv");
    private static final File TARGET_FILE_WITH_MISSING_TARGETS = new File(TEST_SUB_DIR, "agilent_trunc_with_some_missing_targets.tsv");

    private static final SexGenotypeDataCollection TEST_SEX_GENOTYPE_DATA_COLLECTION;

    static {
        try {
            TEST_SEX_GENOTYPE_DATA_COLLECTION = new SexGenotypeDataCollection(TEST_SEX_GENOTYPE_FILE);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the test sex genotype table file: " +
                    TEST_SEX_GENOTYPE_FILE.getAbsolutePath());
        }
    }

    @Test
    public void testGoodRunWithoutInputTargetList() {
        final File sexGenotypeOutputFile = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "-" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypeOutputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        try {
            final SexGenotypeDataCollection inferredSexGenotypeDataCollection =
                    new SexGenotypeDataCollection(sexGenotypeOutputFile);
            Assert.assertEquals(inferredSexGenotypeDataCollection, TEST_SEX_GENOTYPE_DATA_COLLECTION);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
    }

    @Test
    public void testGoodRunWithInputTargetList() {
        final File sexGenotypeOutputFile = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "-" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypeOutputFile.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME, TARGET_FILE_WITH_MISSING_TARGETS.getAbsolutePath()
        };
        runCommandLine(arguments);
        try {
            final SexGenotypeDataCollection inferredSexGenotypeDataCollection =
                    new SexGenotypeDataCollection(sexGenotypeOutputFile);
            Assert.assertEquals(inferredSexGenotypeDataCollection, TEST_SEX_GENOTYPE_DATA_COLLECTION);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
    }

    @Test(dataProvider = "badMappingErrorDataProvider", expectedExceptions = IllegalArgumentException.class)
    public void testBadRunMappingError(final double mappingError) {
        final File sexGenotypeOutputFile = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "-" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypeOutputFile.getAbsolutePath(),
                "-" + TargetCoverageSexGenotyper.BASELINE_MAPPING_ERROR_PROBABILITY_LONG_NAME, String.valueOf(mappingError)
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testBadRunMissingReadCounts() {
        final File sexGenotypeOutputFile = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "-" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypeOutputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testBadRunMissingAnnotations() {
        final File sexGenotypeOutputFile = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypeOutputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testBadRunMissingOutput() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "-" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }

    @DataProvider(name = "badMappingErrorDataProvider")
    public Object[][] badMappingErrorDataProvider() {
        return new Object[][] {{0}, {1}, {-1}, {1.1}};
    }

}
