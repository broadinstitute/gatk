package org.broadinstitute.hellbender.tools.archive;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionIndicator;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionSimulatedData;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountWithPhasePosteriors;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountWithPhasePosteriorsCollection;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Tests for {@link CalculatePulldownPhasePosteriors}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class CalculatePulldownPhasePosteriorsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";
    private static final File TUMOR_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR, "snps-for-allelic-integration.tsv");
    private static final File ACNV_SEGMENT_FILE = new File(TEST_SUB_DIR, "acnv-segments-from-allelic-integration.seg");
    private static final File AF_PARAMS_FILE = new File(TEST_SUB_DIR, "af-params-from-allelic-integration.af.param");
    private static final File ALLELIC_PON_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-normal.tsv");
    private static final String OUTPUT_FILE_NAME = "phase.tsv";

    private static final double FRACTION_OF_INDICATORS_CORRECT_THRESHOLD = 0.95;

    /**
     * Uses {@link AlleleFractionSimulatedData} to test recovery of phase indicators.  Phase with highest posterior
     * probability is compared to the true phase; we require that
     * {@link CalculatePulldownPhasePosteriorsIntegrationTest#FRACTION_OF_INDICATORS_CORRECT_THRESHOLD} of the
     * indicators are recovered correctly.
     */
    @Test
    public void testCalculatePhasePosteriors() {
        final double averageHetsPerSegment = 100;
        final int numSegments = 100;
        final int averageDepth = 100;
        final double biasMean = 1.1;
        final double biasVariance = 0.01;
        final double outlierProbability = 0.02;
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(
                averageHetsPerSegment, numSegments, averageDepth, biasMean, biasVariance, outlierProbability);
        final SegmentedGenome segmentedGenome = simulatedData.getSegmentedGenome();
        final AlleleFractionState trueState = simulatedData.getTrueState();
        final AlleleFractionSimulatedData.PhaseIndicators truePhases = simulatedData.getTruePhases();

        final AllelicCountCollection counts = new AllelicCountCollection();
        segmentedGenome.getGenome().getSNPs().targets().stream().forEach(counts::add);  //note that chromosomes are in lexicographical order
        final AllelicCountWithPhasePosteriorsCollection countsWithPhasePosteriors =
                CalculatePulldownPhasePosteriors.calculatePhasePosteriors(counts, segmentedGenome.getSegments(), trueState, AllelicPanelOfNormals.EMPTY_PON);

        int numIndicatorsCorrect = 0;
        final Iterator<AlleleFractionIndicator> phaseIterator = truePhases.iterator();  //order is ALT_MINOR, REF_MINOR, OUTLIER
        final Iterator<AllelicCountWithPhasePosteriors> countWithPhasePosteriorsIterator = countsWithPhasePosteriors.getCounts().iterator();
        while(phaseIterator.hasNext() && countWithPhasePosteriorsIterator.hasNext()) {
            final AlleleFractionIndicator truePhase = phaseIterator.next();
            final AllelicCountWithPhasePosteriors countWithPhasePosteriors = countWithPhasePosteriorsIterator.next();
            final List<Double> phaseProbabilities = Arrays.asList(
                    countWithPhasePosteriors.getAltMinorProb(),
                    countWithPhasePosteriors.getRefMinorProb(),
                    countWithPhasePosteriors.getOutlierProb());
            final int indexOfMaxProbPhase = phaseProbabilities.indexOf(Collections.max(phaseProbabilities));
            final AlleleFractionIndicator maxProbPhase = AlleleFractionIndicator.values()[indexOfMaxProbPhase];
            if (maxProbPhase.equals(truePhase)) {
                numIndicatorsCorrect++;
            }
        }
        final double fractionOfIndicatorsCorrect = (double) numIndicatorsCorrect / countsWithPhasePosteriors.getCounts().size();
        Assert.assertTrue(fractionOfIndicatorsCorrect >= FRACTION_OF_INDICATORS_CORRECT_THRESHOLD);
    }

    @Test
    public void testRunWithoutAllelicPoN() {
        final File tempDir = createTempDir("phase-posteriors-integration-without-pon");
        final String outputFileName = tempDir.getAbsolutePath() + "/" + OUTPUT_FILE_NAME;

        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.AF_PARAMETER_FILE_LONG_NAME, AF_PARAMS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputFileName,
                "--verbosity", "INFO",
        };
        testRun(arguments, outputFileName);
    }

    @Test
    public void testRunWithAllelicPoN() {
        final File tempDir = createTempDir("phase-posteriors-integration-with-pon");
        final String outputFileName = tempDir.getAbsolutePath() + "/" + OUTPUT_FILE_NAME;

        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, ACNV_SEGMENT_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.AF_PARAMETER_FILE_LONG_NAME, AF_PARAMS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME, ALLELIC_PON_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputFileName,
                "--verbosity", "INFO",
        };
        testRun(arguments, outputFileName);
    }

    private void testRun(final String[] arguments, final String outputFileName) {
        runCommandLine(arguments);

        //only check that file is created and can be read as an AllelicCountWithPhasePosteriorsCollection,
        //do not check for correctness of results (which is tested by testCalculatePhasePosteriors)
        final File outputFile = new File(outputFileName);
        Assert.assertTrue(outputFile.isFile(), outputFile.getAbsolutePath() + " is not a file.");
        Assert.assertTrue(outputFile.length() > 0);
        try {
            //check that file has:
            //  - at least two lines and all either start with "#" or contain at least one "\t"
            //  - at least two lines with tab (column names + 1 SNP)
            final List<String> outputLines = FileUtils.readLines(outputFile);
            Assert.assertTrue(outputLines.size() >= 2);
            Assert.assertEquals(outputLines.stream().filter(l -> l.contains("\t") || l.startsWith("#")).count(), outputLines.size());
            Assert.assertTrue(outputLines.stream().filter(l -> l.split("\t").length > 2 && !l.startsWith("#")).count() > 2,
                    "File: " + outputFile + " does not seem to have at least one SNP and a header.");

            final AllelicCountWithPhasePosteriorsCollection counts = new AllelicCountWithPhasePosteriorsCollection(outputFile);
        } catch (final Exception e) {
            Assert.fail("Could not create AllelicCountWithPhasePosteriorsCollection from file: " + outputFile, e);
        }
    }
}