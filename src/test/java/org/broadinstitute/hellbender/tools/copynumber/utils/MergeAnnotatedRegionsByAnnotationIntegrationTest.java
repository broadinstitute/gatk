package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class MergeAnnotatedRegionsByAnnotationIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "/copynumber/utils/";
    private static final String SIMPLE_TEST_FILE = TEST_SUB_DIR + "merge-annotated-regions-by-annotation-legacy-test.seg";
    private static final String REF = hg19MiniReference;

    /**
     * Very simple test, since this tool basically calls out to a class that does all the heavy lifting.  Most testing
     *  is done in
     * {@link org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalUtilsUnitTest}.
     */
    @Test
    public void testLegacySegFile() {
        final File outputFile = IOUtils.createTempFile("mergeannotatedregionsbyannotation", ".seg");
        final String annotationToMatch = "Segment_Mean";
        final List<String> arguments = new ArrayList<>();
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(SIMPLE_TEST_FILE);
        arguments.add("--" + MergeAnnotatedRegionsByAnnotation.ANNOTATIONS_TO_MATCH);
        arguments.add(annotationToMatch);
        arguments.add("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(collection.getRecords().size(), 2);

        final List<Double> gtSegMeans = Arrays.asList(0.037099, 0.001748);
        final List<Double> segMeans = collection.getRecords().stream().map(r -> r.getAnnotationValue(annotationToMatch))
                .map(Double::parseDouble)
                .collect(Collectors.toList());
        Assert.assertEquals(segMeans, gtSegMeans);
    }
}
