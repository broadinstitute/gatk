package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RevertSamSparkUnitTest extends CommandLineProgramTest {

    private final File basicSamToRevert = getTestFile("revert_sam_basic.sam");

    private final File validOutputMap = getTestFile("revert_sam_valid_output_map.txt");
    private final File nonExistentOutputMap = getTestFile("revert_sam_does_not_exist.txt");
    private final File badHeaderOutputMap = getTestFile("revert_sam_bad_header_output_map.txt");

    @Test
    public static void testFilePathsWithoutMapFile() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");

        final Map<String, Path> outputMap = RevertSamSpark.getOutputMap(null, new File("/foo/bar").getAbsolutePath(), ".bam", Arrays.asList(rg1, rg2), true);
        Assert.assertEquals(outputMap.get("rg1"), IOUtils.getPath(new File("/foo/bar/rg1.bam").getAbsolutePath()));
        Assert.assertEquals(outputMap.get("rg2"), IOUtils.getPath(new File("/foo/bar/rg2.bam").getAbsolutePath()));
    }

    @Test
    public void testValidateOutputParamsByReadGroupMapValid() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsByReadGroup(null, validOutputMap.getAbsolutePath());
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsByReadGroupMissingMap() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsByReadGroup(null, nonExistentOutputMap.getAbsolutePath());
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Cannot read"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupBadHeaderMap() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsByReadGroup(null, badHeaderOutputMap.getAbsolutePath());
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Invalid header"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupNoMapOrDir() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsByReadGroup(null, null);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Must provide either"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupDirValid() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsByReadGroup(createTempDir("testValidateOutputParamsNotByReadGroupValid").getAbsolutePath(), null);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupValid() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsNotByReadGroup(createTempFile("testValidateOutputParamsNotByReadGroupValid","").getAbsolutePath(), null);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupNoOutput() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsNotByReadGroup(null, null);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("output is required"), true);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupMap() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsNotByReadGroup(null, validOutputMap.getAbsolutePath());
        Assert.assertEquals(errors.size(), 2);
        Assert.assertEquals(errors.get(0).contains("Cannot provide outputMap"), true);
        Assert.assertEquals(errors.get(1).contains("output is required"), true);
    }

    @Test
    public static void testGetDefaultExtension() {
        Assert.assertEquals(RevertSamSpark.getDefaultExtension(new GATKPath("this.is.a.sam"), RevertSamSpark.FileType.dynamic), ".sam");
        Assert.assertEquals(RevertSamSpark.getDefaultExtension(new GATKPath("this.is.a.cram"), RevertSamSpark.FileType.dynamic), ".cram"); //TODO https://github.com/broadinstitute/gatk/issues/5559
        Assert.assertEquals(RevertSamSpark.getDefaultExtension(new GATKPath("this.is.a.bam"), RevertSamSpark.FileType.dynamic), ".bam");
        Assert.assertEquals(RevertSamSpark.getDefaultExtension(new GATKPath("foo"), RevertSamSpark.FileType.dynamic), ".bam");
    }

    @Test
    public static void testValidateOutputParamsNotByReadGroupDir() throws IOException {
        final List<String> errors = RevertSamSpark.validateOutputParamsNotByReadGroup(createTempDir("testValidateOutputParamsNotByReadGroupDir").getAbsolutePath(), null);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("should not be a directory"), true);
    }

    @Test
    public void testAssertAllReadGroupsMappedSuccess() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");

        final Map<String, Path> outputMap = new HashMap<>();
        outputMap.put("rg1", IOUtils.getPath(new File("/foo/bar/rg1.bam").getAbsolutePath()));
        outputMap.put("rg2", IOUtils.getPath(new File("/foo/bar/rg2.bam").getAbsolutePath()));
        RevertSamSpark.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1, rg2));
        RevertSamSpark.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1));
        RevertSamSpark.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg2));
    }

    @Test(expectedExceptions = {GATKException.class})
    public void testAssertAllReadGroupsMappedFailure() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");
        final SAMReadGroupRecord rg3 = new SAMReadGroupRecord("rg3");

        final Map<String, Path> outputMap = new HashMap<>();
        outputMap.put("rg1", IOUtils.getPath(new File("/foo/bar/rg1.bam").getAbsolutePath()));
        outputMap.put("rg2", IOUtils.getPath(new File("/foo/bar/rg2.bam").getAbsolutePath()));
        RevertSamSpark.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1, rg2, rg3));
    }
}
