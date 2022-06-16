package org.broadinstitute.hellbender.tools.dragstr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CalibrateDragstrModelIntegrationTest extends CommandLineProgramTest {

    private static final File HUMAN_B37_20_21_TBL_GZ =
            new File("src/test/resources/org/broadinstitute/hellbender/tools/dragstr/human-b37-decimated-20-21.strtbl");
    private static final File EXPECTED_PARAMS_FILE =
            new File("src/test/resources/org/broadinstitute/hellbender/tools/dragstr/expected-ceutrio-hiseq-wgs-b37-20-21.drgstr");

    @DataProvider
    public static Object[][] threadCountAndBamCramData() {
        return new Object[][]{ {1, false}, {1, true}, {3, true}, {3, false}, {4, false}, {4, true}};
    }

    //TODO Ideally the bam/cram would be large enough so that forcing the estimation is not necessary.
    //TODO we should cosider not only human but also other species sample that are nore compect and
    //TODO with larger propr
    @Test(dataProvider = "threadCountAndBamCramData")
    public void testDragstrModelInferenceForcingEstimation(final int threads, final boolean useCram) throws IOException {
        final File output = File.createTempFile("test-output", ".dragstr");
        output.deleteOnExit();

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(StandardArgumentDefinitions.REFERENCE_SHORT_NAME, b37_reference_20_21);
        args.add(CalibrateDragstrModel.STR_TABLE_PATH_FULL_NAME, HUMAN_B37_20_21_TBL_GZ);
        args.addFlag(CalibrateDragstrModel.PARALLEL_FULL_NAME);
        if (threads > 1) args.add(CalibrateDragstrModel.THREADS_FULL_NAME, threads);
        if (threads == 0) args.addFlag(CalibrateDragstrModel.PARALLEL_FULL_NAME);
        args.add(CalibrateDragstrModel.SHARD_SIZE_FULL_NAME,  CalibrateDragstrModel.DEFAULT_SHARD_SIZE);
        args.add(CalibrateDragstrModel.DOWN_SAMPLE_SIZE_FULL_NAME, CalibrateDragstrModel.DEFAULT_DOWN_SAMPLE_SIZE);
        args.addInput(useCram ? NA12878_20_21_WGS_cram : NA12878_20_21_WGS_bam);
        args.addFlag(CalibrateDragstrModel.FORCE_ESTIMATION_FULL_NAME);
        args.addOutput(output);
        final Object res = runCommandLine(args);

        final DragstrParams actual = DragstrParamUtils.parse(new GATKPath(output.toString()));
        final DragstrParams expected = DragstrParamUtils.parse(new GATKPath(EXPECTED_PARAMS_FILE.toString()));
        assertEquals(actual, expected);
    }

    @Test(dataProvider = "threadCountAndBamCramData")
    public void testDragstrModelInferenceFailingOverToDefaults(final int threads, final boolean useCram) throws IOException {
        final File output = File.createTempFile("test-output", ".drgstr");
        output.deleteOnExit();

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(StandardArgumentDefinitions.REFERENCE_SHORT_NAME, b37_reference_20_21);
        args.add(CalibrateDragstrModel.STR_TABLE_PATH_FULL_NAME, HUMAN_B37_20_21_TBL_GZ);
        args.addFlag(CalibrateDragstrModel.PARALLEL_FULL_NAME);
        if (threads > 1) args.add(CalibrateDragstrModel.THREADS_FULL_NAME, threads);
        if (threads == 0) args.addFlag(CalibrateDragstrModel.PARALLEL_FULL_NAME);
        args.add(CalibrateDragstrModel.SHARD_SIZE_FULL_NAME,  CalibrateDragstrModel.DEFAULT_SHARD_SIZE);
        args.add(CalibrateDragstrModel.DOWN_SAMPLE_SIZE_FULL_NAME, CalibrateDragstrModel.DEFAULT_DOWN_SAMPLE_SIZE);
        args.addInput(useCram ? NA12878_20_21_WGS_cram : NA12878_20_21_WGS_bam);
        args.addOutput(output);
        runCommandLine(args);

        final DragstrParams actual = DragstrParamUtils.parse(new GATKPath(output.toString()));
        final DragstrParams expected = DragstrParams.DEFAULT;
        assertEquals(actual, expected);
    }

    private void assertEquals(final DragstrParams actual, final DragstrParams expected) {
        Assert.assertEquals(actual.maximumPeriod(), expected.maximumPeriod());
        Assert.assertEquals(actual.maximumRepeats(), expected.maximumRepeats());
        for (int i = 1; i < actual.maximumPeriod(); i++) {
            for (int j = 1; j < actual.maximumRepeats(); j++) {
                Assert.assertEquals(actual.gop(i, j), expected.gop(i, j), 0.001, "gop (" + i + " , " + j + ") differ");
                Assert.assertEquals(actual.gcp(i, j), expected.gcp(i, j), 0.001, "gcp (" + i + " , " + j + ") differ");
                Assert.assertEquals(actual.api(i, j), expected.api(i, j), 0.001, "api (" + i + " , " + j + ") differ");
            }
        }
    }
}