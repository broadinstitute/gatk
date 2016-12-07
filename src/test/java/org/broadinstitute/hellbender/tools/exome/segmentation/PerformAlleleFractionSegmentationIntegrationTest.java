package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * Created by davidben on 5/23/16.
 */
public final class PerformAlleleFractionSegmentationIntegrationTest extends CommandLineProgramTest {
    private static final String TOOLS_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";
    private static final File ALLELIC_COUNTS_FILE = new File(TOOLS_TEST_DIRECTORY, "snps-for-allelic-integration.tsv");

    @Test
    public void testCommandLine() {
        final File snpFile = ALLELIC_COUNTS_FILE;
        final File outputSegmentFile = createTempFile("segments", ".seg");
        final int initialNumStates = 10;
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, snpFile.getAbsolutePath(),
                "-" + PerformAlleleFractionSegmentation.INITIAL_NUM_STATES_SHORT_NAME, Integer.toString(initialNumStates),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, outputSegmentFile.getAbsolutePath(),
        };
        runCommandLine(arguments);

        final List<ModeledSegment> segments = SegmentUtils.readModeledSegmentsFromSegmentFile(outputSegmentFile);
    }
}