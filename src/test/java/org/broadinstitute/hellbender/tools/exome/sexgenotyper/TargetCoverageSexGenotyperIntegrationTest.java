package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Set;

/**
 * Integration tests for {@link TargetCoverageSexGenotyper}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class TargetCoverageSexGenotyperIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/sexgenotyper/";
    private static final File TEST_CONTIG_PLOIDY_ANNOTS_FILE = new File(TEST_SUB_DIR, "contig_annots.tsv");
    private static final File TEST_RCC_FILE = new File(TEST_SUB_DIR, "ice_trunc_rcc.tsv");
    private static final File TEST_SEX_GENOTYPE_FILE = new File(TEST_SUB_DIR, "ice_trunc_sex_genotypes.tsv");
    private static final File TARGET_FILE_WITH_MISSING_TARGETS = new File(TEST_SUB_DIR, "ice_trunc_targets_some_missing.tsv");
    private static final File TARGET_FILE_ALL_TARGETS_WITH_BAIT_COUNTS = new File(TEST_SUB_DIR, "ice_trunc_targets_with_bait_counts.tsv");

    private static final double EPS = 1e-12;

    private File sexGenotypesWithMissingTargets;
    private File sexGenotypesWithExcludedIntervals;
    private File sexGenotypesWithoutBaitCounts;
    private File sexGenotypesWithBaitCounts;

    private static final SexGenotypeDataCollection TEST_SEX_GENOTYPE_DATA_COLLECTION;

    static {
        try {
            TEST_SEX_GENOTYPE_DATA_COLLECTION = new SexGenotypeDataCollection(TEST_SEX_GENOTYPE_FILE);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the test sex genotype table file: " +
                    TEST_SEX_GENOTYPE_FILE.getAbsolutePath());
        }
    }

    /**
     * Basic run
     */
    @Test
    public void testSuccessfulRunWithoutInputTargetList() {
        sexGenotypesWithoutBaitCounts = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypesWithoutBaitCounts.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
        try {
            final SexGenotypeDataCollection inferredSexGenotypeDataCollection =
                    new SexGenotypeDataCollection(sexGenotypesWithoutBaitCounts);
            Assert.assertEquals(inferredSexGenotypeDataCollection, TEST_SEX_GENOTYPE_DATA_COLLECTION);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
    }

    /**
     * A user can provide a genotyping target list. The target list can include bait counts to correct for
     * the associated coverage bias.
     */
    @Test
    public void testSuccessfulRunWithBaitCounts() {
        sexGenotypesWithBaitCounts = createTempFile("sex-genotype-with-bait-counts", ".tsv");
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypesWithBaitCounts.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME, TARGET_FILE_ALL_TARGETS_WITH_BAIT_COUNTS.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
        try {
            final SexGenotypeDataCollection inferredSexGenotypeDataCollection =
                    new SexGenotypeDataCollection(sexGenotypesWithBaitCounts);
            Assert.assertEquals(inferredSexGenotypeDataCollection, TEST_SEX_GENOTYPE_DATA_COLLECTION);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
    }

    /**
     * Certain targets can be excluded form the read count table by providing a custom target list to the tool. If only
     * a fraction of targets are excluded, inferred sex genotypes are not affected (though the likelihoods will change).
     *
     * Here, the input target list excludes chr1 as well as X:14708580-17152282.
     */
    @Test
    public void testSuccessfulRunWithMissingTargetsWithoutBaitsCounts() {
        sexGenotypesWithMissingTargets = createTempFile("sex-genotype-missing-targets", ".tsv");
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypesWithMissingTargets.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME, TARGET_FILE_WITH_MISSING_TARGETS.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
        try {
            final SexGenotypeDataCollection inferredSexGenotypeDataCollection =
                    new SexGenotypeDataCollection(sexGenotypesWithMissingTargets);
            Assert.assertEquals(inferredSexGenotypeDataCollection, TEST_SEX_GENOTYPE_DATA_COLLECTION);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
    }

    /**
     * Targets can also be excluded via command line. If only a fraction of targets are excluded,
     * inferred sex genotypes are not affected (though the likelihoods will change).
     */
    @Test
    public void testSuccessfulRunWithExcludedIntervals() {
        sexGenotypesWithExcludedIntervals = createTempFile("sex-genotype-excluded-intervals", ".tsv");
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypesWithExcludedIntervals.getAbsolutePath(),
                "--excludeIntervals", "1",
                "--excludeIntervals", "X:14708580-17152282",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
        try {
            final SexGenotypeDataCollection inferredSexGenotypeDataCollection =
                    new SexGenotypeDataCollection(sexGenotypesWithExcludedIntervals);
            Assert.assertEquals(inferredSexGenotypeDataCollection, TEST_SEX_GENOTYPE_DATA_COLLECTION);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
    }

    @Test(dependsOnMethods = {"testSuccessfulRunWithMissingTargetsWithoutBaitsCounts", "testSuccessfulRunWithExcludedIntervals"})
    public void testIntervalExclusionSameBehaviorAsTrimmingTargetList() {
        Assert.assertNotNull(sexGenotypesWithMissingTargets);
        Assert.assertNotNull(sexGenotypesWithExcludedIntervals);
        final SexGenotypeDataCollection sexGenotypesCollectionWithMissingTargets;
        final SexGenotypeDataCollection sexGenotypesCollectionWithExcludedIntervals;
        try {
            sexGenotypesCollectionWithMissingTargets = new SexGenotypeDataCollection(sexGenotypesWithMissingTargets);
            sexGenotypesCollectionWithExcludedIntervals = new SexGenotypeDataCollection(sexGenotypesWithExcludedIntervals);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
        Assert.assertTrue(checkSexGenotypeDataCollectionCompleteEquality(sexGenotypesCollectionWithMissingTargets,
                sexGenotypesCollectionWithExcludedIntervals));
    }

    /**
     * A target list that includes BAIT_COUNT annotations result in different (and more reliable) sex genotype likelihoods
     */
    @Test(dependsOnMethods = {"testSuccessfulRunWithoutInputTargetList", "testSuccessfulRunWithBaitCounts"})
    public void testBaitCountsChangeGenotypeLikelihood() {
        Assert.assertNotNull(sexGenotypesWithoutBaitCounts);
        Assert.assertNotNull(sexGenotypesWithBaitCounts);
        final SexGenotypeDataCollection sexGenotypesCollectionWithoutBaitCounts;
        final SexGenotypeDataCollection sexGenotypesCollectionWithBaitCounts;
        try {
            sexGenotypesCollectionWithoutBaitCounts = new SexGenotypeDataCollection(sexGenotypesWithoutBaitCounts);
            sexGenotypesCollectionWithBaitCounts = new SexGenotypeDataCollection(sexGenotypesWithBaitCounts);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
        Assert.assertFalse(checkSexGenotypeDataCollectionCompleteEquality(sexGenotypesCollectionWithoutBaitCounts,
                sexGenotypesCollectionWithBaitCounts));
    }

    /**
     * It is likely to be able to determine the correct sex genotype even without ChrY coverage (just based on ChrX ploidy
     * difference between XX and XY). This, however, requires a more reliable modeling of coverage biases. Taking
     * into account bait count bias is a necessary ingredient (see the next test).
     */
    @Test
    public void testSuccessfulRunWithBaitCountsWithoutChrY() {
        final File sexGenotypesWithBaitCountsWithoutChrY = createTempFile("sex-genotype-with-bait-counts-without-chrY", ".tsv");
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypesWithBaitCountsWithoutChrY.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME, TARGET_FILE_ALL_TARGETS_WITH_BAIT_COUNTS.getAbsolutePath(),
                "--excludeIntervals", "Y",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
        try {
            final SexGenotypeDataCollection inferredSexGenotypeDataCollection =
                    new SexGenotypeDataCollection(sexGenotypesWithBaitCountsWithoutChrY);
            Assert.assertEquals(inferredSexGenotypeDataCollection, TEST_SEX_GENOTYPE_DATA_COLLECTION);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
    }

    /**
     * Without ChrY coverage and correcting for bait count bias, sex genotyping can be unreliable. The model is biased
     * toward higher X-ploidy sex genotypes (here, XX) if the X contig has a few targets with anomalously large number
     * of baits.
     */
    @Test
    public void testSuccessfulRunWithoutBaitCountsWithoutChrY() {
        final File sexGenotypesWithoutBaitCountsWithoutChrY = createTempFile("sex-genotype-without-bait-counts-without-chrY", ".tsv");
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypesWithoutBaitCountsWithoutChrY.getAbsolutePath(),
                "--excludeIntervals", "Y",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
        try {
            final SexGenotypeDataCollection inferredSexGenotypeDataCollection =
                    new SexGenotypeDataCollection(sexGenotypesWithoutBaitCountsWithoutChrY);
            /* we assert inequality */
            Assert.assertNotEquals(inferredSexGenotypeDataCollection, TEST_SEX_GENOTYPE_DATA_COLLECTION);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the output sex genotype table");
        }
    }

    /**
     * Asserts that two sex genotype data collections are completely equal (sample names, sex genotypes, likelihoods)
     */
    private static boolean checkSexGenotypeDataCollectionCompleteEquality(final SexGenotypeDataCollection first,
                                                                          final SexGenotypeDataCollection second) {
        if (first == null || second == null) {
            return false;
        }
        if (!first.getSampleNames().equals(second.getSampleNames())) {
            return false;
        }
        final Set<String> sampleNames = first.getSampleNames();
        for (final String sampleName : sampleNames) {
            final SexGenotypeData firstSexGenotypeData = first.getSampleSexGenotypeData(sampleName);
            final SexGenotypeData secondSexGenotypeData = second.getSampleSexGenotypeData(sampleName);
            if (!firstSexGenotypeData.getSexGenotypesSet().equals(secondSexGenotypeData.getSexGenotypesSet())) {
                return false;
            }
            final Set<String> sexGenotypes = firstSexGenotypeData.getSexGenotypesSet();
            for (final String sexGenotype : sexGenotypes) {
                if (FastMath.abs(firstSexGenotypeData.getLogLikelihoodPerGenotype(sexGenotype) -
                        secondSexGenotypeData.getLogLikelihoodPerGenotype(sexGenotype)) > EPS) {
                    return false;
                }
            }
        }
        return true;
    }

    @Test(dataProvider = "badMappingErrorDataProvider", expectedExceptions = IllegalArgumentException.class)
    public void testBadRunMappingError(final double mappingError) {
        final File sexGenotypeOutputFile = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypeOutputFile.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.BASELINE_MAPPING_ERROR_PROBABILITY_LONG_NAME, String.valueOf(mappingError),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testBadRunMissingReadCounts() {
        final File sexGenotypeOutputFile = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypeOutputFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testBadRunMissingAnnotations() {
        final File sexGenotypeOutputFile = createTempFile("sex-genotype-test", ".tsv");
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, sexGenotypeOutputFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testBadRunMissingOutput() {
        final String[] arguments = {
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, TEST_RCC_FILE.getAbsolutePath(),
                "--" + TargetCoverageSexGenotyper.INPUT_CONTIG_ANNOTATIONS_LONG_NAME, TEST_CONTIG_PLOIDY_ANNOTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @DataProvider(name = "badMappingErrorDataProvider")
    public Object[][] badMappingErrorDataProvider() {
        return new Object[][] {{0}, {1}, {-1}, {1.1}};
    }
}
