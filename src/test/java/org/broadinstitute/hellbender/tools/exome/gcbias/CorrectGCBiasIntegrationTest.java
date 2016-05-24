package org.broadinstitute.hellbender.tools.exome.gcbias;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.fakedata.GCBiasSimulatedData;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class CorrectGCBiasIntegrationTest extends CommandLineProgramTest {
    // Test meta-parameters:
    private static final int NUM_TARGETS = 10000;
    private static final int NUM_SAMPLES = 10;

    private static final File TARGETS_FILE = createTempFile("targets", ".tsv");
    private static final File INPUT_COUNTS_FILE = createTempFile("input", ".tsv");
    private static final File OUTPUT_COUNTS_FILE = createTempFile("output", ".tsv");

    double[] gcContentByTarget;
    ReadCountCollection inputCounts;


    @BeforeClass
    public void createInputFiles() throws IOException {
        final Pair<ReadCountCollection, double[]> data = GCBiasSimulatedData.simulatedData(NUM_TARGETS, NUM_SAMPLES);
        inputCounts = data.getLeft();
        gcContentByTarget = data.getRight();
        GCBiasSimulatedData.makeGCBiasInputFiles(data, INPUT_COUNTS_FILE, TARGETS_FILE);
    }

    // test that results match expected behavior of the backing class
    @Test
    public void testGCCorrection() throws IOException {
        final List<String> arguments = new ArrayList<>();
        arguments.addAll(Arrays.asList(
                "-" + CorrectGCBias.INPUT_READ_COUNTS_FILE_SHORT_NAME, INPUT_COUNTS_FILE.getAbsolutePath(),
                "-" + CorrectGCBias.OUTPUT_READ_COUNTS_FILE_SHORT_NAME, OUTPUT_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TARGETS_FILE.getAbsolutePath()));
        runCommandLine(arguments);

        final ReadCountCollection outputCounts = ReadCountCollectionUtils.parse(OUTPUT_COUNTS_FILE);
        final ReadCountCollection expectedOutputCounts = GCCorrector.correctCoverage(inputCounts, gcContentByTarget);
        Assert.assertEquals(outputCounts.columnNames(), inputCounts.columnNames());
        Assert.assertEquals(outputCounts.counts().subtract(expectedOutputCounts.counts()).getNorm(), 0, 1e-10);
    }

    @AfterClass
    public void disposeFile() throws IOException {
        TARGETS_FILE.delete();
        INPUT_COUNTS_FILE.delete();
        OUTPUT_COUNTS_FILE.delete();
    }
}