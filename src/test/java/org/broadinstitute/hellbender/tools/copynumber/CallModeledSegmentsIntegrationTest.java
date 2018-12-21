package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledModeledSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledModeledSegment;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.io.File;

public class CallModeledSegmentsIntegrationTest extends CommandLineProgramTest {
    // Simulated samples
    private static final File SIMULATED_DATA_DIR = new File(GATKBaseTest.publicTestDir,"org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/");
    private static final File OUTPUT_DIR = createTempDir("out");
    private static final File INPUT_SIMULATED_NORMAL_DATA = new File(SIMULATED_DATA_DIR, "normal_sample_data.seg");
    private static final File INPUT_SIMULATED_NORMAL_TRUTH = new File(SIMULATED_DATA_DIR, "normal_sample_truth.seg");
    private static final File INPUT_SIMULATED_NORMAL_NOISY_DATA = new File(SIMULATED_DATA_DIR, "normal_noisy_sample_data.seg");
    private static final File INPUT_SIMULATED_NORMAL_NOISY_TRUTH = new File(SIMULATED_DATA_DIR, "normal_noisy_sample_truth.seg");
    private static final File INPUT_SIMULATED_SOMATIC_40P_PURITY_DATA = new File(SIMULATED_DATA_DIR, "cancer_0p40_purity_data.seg");
    private static final File INPUT_SIMULATED_SOMATIC_40P_PURITY_TRUTH = new File(SIMULATED_DATA_DIR, "cancer_0p40_purity_truth.seg");
    private static final File INPUT_SIMULATED_SOMATIC_60P_PURITY_DATA = new File(SIMULATED_DATA_DIR, "cancer_0p60_purity_data.seg");
    private static final File INPUT_SIMULATED_SOMATIC_60P_PURITY_TRUTH = new File(SIMULATED_DATA_DIR, "cancer_0p60_purity_truth.seg");
    private static final File INPUT_SIMULATED_SOMATIC_100P_PURITY_DATA = new File(SIMULATED_DATA_DIR, "cancer_1p00_purity_data.seg");
    private static final File INPUT_SIMULATED_SOMATIC_100P_PURITY_TRUTH = new File(SIMULATED_DATA_DIR, "cancer_1p00_purity_truth.seg");

    @DataProvider
    public Iterator<Object[]> simulatedData() {
        // The data provider contains the inputs of the runTest method
        // inputFile, truthFile, loadCopyRatio, loadAlleleFraction, interactiveMode

        // TESTING: two inputs, for which the tests didn't pass, were commented out below.
        // We shall restore this later.
        List<Object[]> result = new LinkedList<>();
        result.add(new Object[] {INPUT_SIMULATED_NORMAL_DATA, INPUT_SIMULATED_NORMAL_TRUTH,
                true, true, true, "normal_simulated"});
        result.add(new Object[] {INPUT_SIMULATED_NORMAL_NOISY_DATA, INPUT_SIMULATED_NORMAL_NOISY_TRUTH,
                true, true, true, "normal_noisy_simulated"});
        result.add(new Object[] {INPUT_SIMULATED_SOMATIC_40P_PURITY_DATA, INPUT_SIMULATED_SOMATIC_40P_PURITY_TRUTH,
                true, true, true, "somatic_100per_cent"});
        result.add(new Object[] {INPUT_SIMULATED_SOMATIC_60P_PURITY_DATA, INPUT_SIMULATED_SOMATIC_60P_PURITY_TRUTH,
                true, true, true, "somatic_100per_cent"});
        result.add(new Object[] {INPUT_SIMULATED_SOMATIC_100P_PURITY_DATA, INPUT_SIMULATED_SOMATIC_100P_PURITY_TRUTH,
                true, true, true, "somatic_100per_cent"});
        return result.iterator();
    }

    @Test(groups={"python"}, dataProvider = "simulatedData")
    public void testSimulatedData(final File inputFile,
                                  final File truthFile,
                                  final boolean loadCopyRatio,
                                  final boolean loadAlleleFraction,
                                  final boolean interactiveRun,
                                  final String outputPrefix) {
        runTest(inputFile, truthFile,
                loadCopyRatio, loadAlleleFraction, interactiveRun,
                outputPrefix);
    }

    private void runTest(final File inputFile,
                         final File truthFile,
                         final boolean loadCopyRatio,
                         final boolean loadAlleleFraction,
                         final boolean interactiveRun,
                         final String outputPrefix) {
        // Method running the tests.

        // Make sure that either the copy ratio or the allele fraction data is loaded
        if (!loadCopyRatio && !loadAlleleFraction) {
            System.out.println("Error: either copy ratio or allele fraction or both need to be loaded.");
            Assert.assertTrue(false);
        }

        // Run clustering algorithm
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputFile.getAbsolutePath())
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .addArgument(CallModeledSegments.OUTPUT_PREFIX_LONG_NAME, outputPrefix)
                .addArgument(CallModeledSegments.LOAD_COPY_RATIO_LONG_NAME, String.valueOf(loadCopyRatio))
                .addArgument(CallModeledSegments.LOAD_ALLELE_FRACTION_LONG_NAME, String.valueOf(loadAlleleFraction))
                .addArgument(CallModeledSegments.INTERACTIVE_RUN_LONG_NAME, String.valueOf(interactiveRun))
                .addArgument(CallModeledSegments.LOG_LONG_NAME, String.valueOf(true))
                .addArgument(CallModeledSegments.N_INFERENCE_ITERATIONS, String.valueOf(40000));

        String outputCallsFilePath = (OUTPUT_DIR.getAbsolutePath() + "/" + outputPrefix
                + CallModeledSegments.OUTPUT_CALLS_SUFFIX_DEFAULT_VALUE);
        String outputImageFilePath = (OUTPUT_DIR.getAbsolutePath() + "/" + outputPrefix
                + CallModeledSegments.OUTPUT_IMAGE_SUFFIX_DEFAULT_VALUE);
        File outputCallsFile = new File(outputCallsFilePath);
        File outputImageFile = new File(outputImageFilePath);
        runCommandLine(argsBuilder);

        // Make sure that output image file and calls file exist
        assertOutputFiles(outputCallsFile, outputImageFile, OUTPUT_DIR);

        // Make sure that the output calls agree with the truth data
        boolean onlyNormalSegments = compareCalledFiles(outputCallsFile, truthFile);

        // Make sure that the interactive images were indeed produced in interactive run
        if (interactiveRun) {
            assertInteractiveFiles(OUTPUT_DIR, outputPrefix, onlyNormalSegments);
        }
    }

    private static void assertOutputFiles(final File outputCallsFile,
                                          final File outputImageFile,
                                          final File outputDir) {
        // Make sure that these files and folders all exist.
        Assert.assertTrue(outputDir.isDirectory());
        Assert.assertTrue(outputCallsFile.isFile());
        Assert.assertTrue(outputImageFile.isFile());
    }

    private static void assertInteractiveFiles(final File outputImageDir,
                                               final String outputImagePrefix,
                                               boolean onlyNormalSegments) {
        // Make sure that all images were generated, including the ones in interactive mode.
        // Input onlyNormalSegments is true iff the algorithm only found a normal segments.
        Assert.assertTrue(outputImageDir.isDirectory());

        // Images generated in interactive mode:
        // Image showing deletions, amplifications and CNLOH.
        final File delAmplPlot = new File(outputImageDir, outputImagePrefix
                + CallModeledSegments.INTERACTIVE_OUTPUT_CLASSIFICATION_IMAGE_SUFFIX_DEFAULT_VALUE);
        // Summary plot showing copy ratio, allele fraction histograms as well as the the segments in
        // (copy ratio, allele fraction plane), indicating which ones are normal
        final File scatterPlot = new File(outputImageDir, outputImagePrefix
                + CallModeledSegments.INTERACTIVE_OUTPUT_SUMMARY_PLOT_SUFFIX_DEFAULT_VALUE);
        // 1D fit to copy ratio data
        final File copyRatioFitPlot = new File(outputImageDir, outputImagePrefix
                + CallModeledSegments.INTERACTIVE_OUTPUT_COPY_RATIO_PLOT_SUFFIX_DEFAULT_VALUE);

        Assert.assertTrue(delAmplPlot.isFile());
        Assert.assertTrue(scatterPlot.isFile());
        Assert.assertTrue(copyRatioFitPlot.isFile());

        if (!onlyNormalSegments) {
            final File alleleFractionPlot = new File(outputImageDir.getAbsolutePath() + "/" + outputImagePrefix
                    + CallModeledSegments.INTERACTIVE_OUTPUT_ALLELE_FRACTION_PLOT_SUFFIX_DEFAULT_VALUE);
            final File copyRatioClusteringPlot = new File(outputImageDir.getAbsolutePath() + "/" + outputImagePrefix
                    + CallModeledSegments.INTERACTIVE_OUTPUT_COPY_RATIO_CLUSTERING_SUFFIX_DEFAULT_VALUE);

            Assert.assertTrue(alleleFractionPlot.isFile());
            Assert.assertTrue(copyRatioClusteringPlot.isFile());
        }
    }

    private static boolean compareCalledFiles(final File outputFile, final File truthFile) {
        // Make sure that the truth file's and the output file's calls disagree at most in errorTolerance fraction of
        // the base pairs, where errorTolerances is specified as:
        double errorTolerance = 0.025;

        CalledModeledSegmentCollection outputData = new CalledModeledSegmentCollection(outputFile);
        CalledModeledSegmentCollection truthData = new CalledModeledSegmentCollection(truthFile);
        List<CalledModeledSegment> outputDataRecords = outputData.getRecords();
        List<CalledModeledSegment> truthDataRecords = truthData.getRecords();

        // Calculating the frequency of disagreeing base pairs between the output and the truth data
        boolean onlyNormaSegments = true;   // True iff the algorithm has found normal segments only
        double totalLength = 0.;
        double segmentLength;
        double lengthDifferent = 0.;
        for (int i=0; i<outputDataRecords.size(); i++) {
            segmentLength = outputDataRecords.get(i).getEnd() - outputDataRecords.get(i).getStart() + 1;
            totalLength += segmentLength;
            if (!(outputDataRecords.get(i).getCallNormal()).equals(truthDataRecords.get(i).getCallNormal())) {
                lengthDifferent += segmentLength;
            }
            if (!outputDataRecords.get(i).getCallNormal().equals("0")) {
                onlyNormaSegments = false;
            }
        }

        if (lengthDifferent / totalLength > errorTolerance) {
            System.out.println("Error: the calls in the truth data and the output are different:");
            printComparison(outputDataRecords, truthDataRecords);
            Assert.assertTrue(false);
        }

        return onlyNormaSegments;
    }

    private static void printComparison(final List<CalledModeledSegment> outputDataRecords,
                                        final List<CalledModeledSegment> truthDataRecords) {
        // Auxiliary method to print the comparison between the output calls file and the truth file.
        // Lines with identical calls are shown in green, whereas the disagreeing ones are shown in red.
        final String ANSI_RESET = "\u001B[0m";
        final String ANSI_RED = "\u001B[31m";
        final String ANSI_GREEN = "\u001B[32m";
        String outputLine;

        System.out.println("CONTIG\tSTART\tEND\tCOPY_RATIO_MEDIAN\tMINOR_ALLELE_FRACTION_MEDIAN\tTRUTH_CALL\tOUTPUT_CALL\tTRUTH_PHRED\tOUTPUT_PHRED");
        for (int i=0; i<outputDataRecords.size(); i++) {
            outputLine = String.valueOf(outputDataRecords.get(i).getContig()) + "\t"
                        + String.valueOf(outputDataRecords.get(i).getStart()) + "\t"
                        + String.valueOf(outputDataRecords.get(i).getEnd()) + "\t"
                        + String.valueOf(Math.pow(2., outputDataRecords.get(i).getLog2CopyRatioSimplePosteriorSummary().getDecile50())) + "\t"
                        + String.valueOf(outputDataRecords.get(i).getMinorAlleleFractionSimplePosteriorSummary().getDecile50()) + "\t"
                        + String.valueOf(outputDataRecords.get(i).getCallNormal()) + "\t"
                        + String.valueOf(truthDataRecords.get(i).getCallNormal()) + "\t"
                        + String.valueOf(outputDataRecords.get(i).getPHREDScoreNormal()) + "\t"
                        + String.valueOf(truthDataRecords.get(i).getPHREDScoreNormal()) + "\t";
            if (outputDataRecords.get(i).getCallNormal().equals(truthDataRecords.get(i).getCallNormal())) {
                System.out.println(ANSI_GREEN + outputLine + ANSI_RESET);
            } else {
                System.out.println(ANSI_RED + outputLine + ANSI_RESET);
            }
        }
    }
}