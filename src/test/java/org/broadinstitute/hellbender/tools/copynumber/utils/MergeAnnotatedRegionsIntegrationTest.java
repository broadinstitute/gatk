package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.ImmutableSortedMap;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MergeAnnotatedRegionsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/utils/";
    private static final String SIMPLE_TEST_FILE = TEST_SUB_DIR + "merge-annotated-regions-simple-test.seg";
    private static final String REF = hg19MiniReference;

    @Test
    public void basicTest() throws IOException {
        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("mergeannotatedregions", ".seg");
        final List<String> arguments = new ArrayList<>();
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(SIMPLE_TEST_FILE);
        arguments.add("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        final AnnotatedIntervalCollection collection =
                AnnotatedIntervalCollection.create(outputFile.toPath(), null);

        Assert.assertEquals(collection.getRecords().size(), 5);

        // Painstakingly made by hand
        Assert.assertEquals(collection.getRecords().get(0), new AnnotatedInterval(new SimpleInterval("1", 100, 105),
                ImmutableSortedMap.of("Num_Probes", "100", "Segment_Mean", "-0.03", "Segment_Call", "0")));
        Assert.assertEquals(collection.getRecords().get(1), new AnnotatedInterval(new SimpleInterval("1", 525, 2305),
                ImmutableSortedMap.of("Num_Probes", "200", "Segment_Mean", "-0.10__-0.76", "Segment_Call", "-__0")));
        Assert.assertEquals(collection.getRecords().get(2), new AnnotatedInterval(new SimpleInterval("1", 2306, 2605),
                ImmutableSortedMap.of("Num_Probes", "200__500", "Segment_Mean", "-0.60", "Segment_Call", "-")));
        Assert.assertEquals(collection.getRecords().get(3), new AnnotatedInterval(new SimpleInterval("2", 525, 1097),
                ImmutableSortedMap.of("Num_Probes", "200", "Segment_Mean", "-0.76", "Segment_Call", "-")));
        Assert.assertEquals(collection.getRecords().get(4), new AnnotatedInterval(new SimpleInterval("2", 1098, 2305),
                ImmutableSortedMap.of("Num_Probes", "200", "Segment_Mean", "-0.10", "Segment_Call", "0")));
    }
}
