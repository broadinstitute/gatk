package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MergeAnnotatedRegionsByAnnotationIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/utils/";
    private static final String SIMPLE_TEST_FILE = TEST_SUB_DIR + "merge-annotated-regions-by-annotation-legacy-test.seg";
    private static final String REF = hg19MiniReference;

    /**
     * Very simple test, since this tool basically calls out to a class that does all the heavy lifting.
     * @throws IOException
     */
    @Test
    public void testLegacySegFile() throws IOException {
        final File outputFile = File.createTempFile("mergeannotatedregionsbyannotation", ".seg");
        final List<String> arguments = new ArrayList<>();
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(SIMPLE_TEST_FILE);
        arguments.add("--" + MergeAnnotatedRegionsByAnnotation.ANNOTATIONS_TO_MATCH);
        arguments.add("Segment_Mean");
        arguments.add("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(collection.getRecords().size(), 2);

        // TODO: Assert that the values are okay
    }

}
