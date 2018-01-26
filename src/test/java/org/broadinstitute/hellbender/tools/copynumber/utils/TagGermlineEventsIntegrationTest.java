package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegionCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class TagGermlineEventsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_RESOURCE_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/utils/";
    public static final String TAG_GERMLINE_NORMAL_SEG_FILE = TEST_RESOURCE_DIR + "tag-germline-with-legacy-header-normal.tsv";
    public static final String TAG_GERMLINE_TUMOR_MATCHED_NORMAL_SEG_FILE = TEST_RESOURCE_DIR + "tag-germline-with-legacy-header-tumor-all-match-normal.tsv";
    public static final String TAG_GERMLINE_TUMOR_ALMOST_MATCHED_NORMAL_SEG_FILE = TEST_RESOURCE_DIR + "tag-germline-with-legacy-header-tumor-some-match-normal.tsv";
    public static final String TAG_GERMLINE_TUMOR_NOT_MATCHED_NORMAL_SEG_FILE = TEST_RESOURCE_DIR + "tag-germline-with-legacy-header-no-match-normal.tsv";
    public static final String TAG_GERMLINE_TUMOR_SPLIT_ALMOST_MATCHED_NORMAL_SEG_FILE = TEST_RESOURCE_DIR + "tag-germline-with-legacy-header-tumor-split-some-match-normal.tsv";
    public static final String TAG_GERMLINE_TUMOR_SPLIT_NO_MATCHED_NORMAL_SEG_FILE = TEST_RESOURCE_DIR + "tag-germline-with-legacy-header-tumor-split-no-match-normal.tsv";

    public static final String REF = hg19_chr1_1M_Reference;

    @Test
    public void testBasic() throws IOException {
        final File outputFile = File.createTempFile("tag_germline_seg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(TAG_GERMLINE_TUMOR_MATCHED_NORMAL_SEG_FILE);
        arguments.add("--" + TagGermlineEvents.MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME);
        arguments.add(TAG_GERMLINE_NORMAL_SEG_FILE);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(outputFile).getRecords();

        // Test that the germline calls are 0, -, 0, -
        Assert.assertEquals(regions.get(0).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(1).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.DELETION.getOutputString());
        Assert.assertEquals(regions.get(2).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(3).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.DELETION.getOutputString());

        final List<SimpleAnnotatedGenomicRegion> regionsInput = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(new File(TAG_GERMLINE_TUMOR_MATCHED_NORMAL_SEG_FILE)).getRecords();
        assertNoRegionChanges(regions, regionsInput);
    }

    @Test
    public void testSlightlyDifferent() throws IOException {
        final File outputFile = File.createTempFile("tag_germline_seg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(TAG_GERMLINE_TUMOR_ALMOST_MATCHED_NORMAL_SEG_FILE);
        arguments.add("--" + TagGermlineEvents.MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME);
        arguments.add(TAG_GERMLINE_NORMAL_SEG_FILE);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(outputFile).getRecords();

        // Test that the germline calls are 0, 0, 0, -
        Assert.assertEquals(regions.get(0).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(1).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(2).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(3).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.DELETION.getOutputString());

        final List<SimpleAnnotatedGenomicRegion> regionsInput = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(new File(TAG_GERMLINE_TUMOR_ALMOST_MATCHED_NORMAL_SEG_FILE)).getRecords();
        assertNoRegionChanges(regions, regionsInput);
    }

    @Test
    public void testSlightlyDifferentAndSplit() throws IOException {
        final File outputFile = File.createTempFile("tag_germline_seg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(TAG_GERMLINE_TUMOR_SPLIT_ALMOST_MATCHED_NORMAL_SEG_FILE);
        arguments.add("--" + TagGermlineEvents.MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME);
        arguments.add(TAG_GERMLINE_NORMAL_SEG_FILE);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(outputFile).getRecords();

        // Test that the germline calls are 0, 0, 0, -, -
        Assert.assertEquals(regions.get(0).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(1).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(2).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(3).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.DELETION.getOutputString());
        Assert.assertEquals(regions.get(4).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.DELETION.getOutputString());

        final List<SimpleAnnotatedGenomicRegion> regionsInput = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(new File(TAG_GERMLINE_TUMOR_SPLIT_ALMOST_MATCHED_NORMAL_SEG_FILE)).getRecords();
        assertNoRegionChanges(regions, regionsInput);
    }

    @Test
    public void testTotallyDifferent() throws IOException {
        final File outputFile = File.createTempFile("tag_germline_seg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(TAG_GERMLINE_TUMOR_NOT_MATCHED_NORMAL_SEG_FILE);
        arguments.add("--" + TagGermlineEvents.MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME);
        arguments.add(TAG_GERMLINE_NORMAL_SEG_FILE);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(outputFile).getRecords();

        // Test that the germline calls are 0, 0, 0, 0
        Assert.assertEquals(regions.get(0).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(1).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(2).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(3).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());

        final List<SimpleAnnotatedGenomicRegion> regionsInput = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(new File(TAG_GERMLINE_TUMOR_NOT_MATCHED_NORMAL_SEG_FILE)).getRecords();
        assertNoRegionChanges(regions, regionsInput);
    }

    // Does not care if output is a superset of input, but the input must be fully and exactly contained in the output.
    private void assertNoRegionChanges(List<SimpleAnnotatedGenomicRegion> regionsOutput, List<SimpleAnnotatedGenomicRegion> regionsInput) {
        // Make sure that no intervals were changed.
        Assert.assertEquals(regionsOutput.stream().map(SimpleAnnotatedGenomicRegion::getInterval).collect(Collectors.toList()),
                regionsInput.stream().map(SimpleAnnotatedGenomicRegion::getInterval).collect(Collectors.toList()));

        // Make sure that no annotation values were changed (just the one added)
        for (int i = 0; i < regionsInput.size(); i++) {
            Assert.assertTrue(regionsOutput.get(i).getAnnotations().entrySet().containsAll(regionsInput.get(i).getAnnotations().entrySet()));
        }
    }

    @Test
    public void testSplit() throws IOException {
        final File outputFile = File.createTempFile("tag_germline_seg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(TAG_GERMLINE_TUMOR_SPLIT_NO_MATCHED_NORMAL_SEG_FILE);
        arguments.add("--" + TagGermlineEvents.MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME);
        arguments.add(TAG_GERMLINE_NORMAL_SEG_FILE);
        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> regions = SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(outputFile).getRecords();

        // Test that the germline calls are 0, 0, 0, 0, 0
        Assert.assertEquals(regions.get(0).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(1).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(2).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(3).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        Assert.assertEquals(regions.get(4).getAnnotationValue(TagGermlineEvents.GERMLINE_TAG_HEADER), CalledCopyRatioSegment.Call.NEUTRAL.getOutputString());
        final List<SimpleAnnotatedGenomicRegion> regionsInput =SimpleAnnotatedGenomicRegionCollection.readAnnotatedRegions(new File(TAG_GERMLINE_TUMOR_SPLIT_NO_MATCHED_NORMAL_SEG_FILE)).getRecords();
        assertNoRegionChanges(regions, regionsInput);
    }
}