package org.broadinstitute.hellbender.utils.samples;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * Tests for SampleUtils
 */
public class SampleUtilsUnitTest extends BaseTest {

    @Override
    public String getTestedClassName() { return "SampleUtils"; }

    @Test
    public void testReadOneSampleFile() {
        Set<File> sampleFiles = new HashSet<>(1);
        sampleFiles.add(new File(getToolTestDataDir(), "samples1.samples"));
        Collection<String> samples = SampleUtils.getSamplesFromFiles(sampleFiles);
        Assert.assertEquals(samples.size(), 2);
        Assert.assertTrue(samples.contains("NA20845"));
        Assert.assertTrue(samples.contains("NA20846"));
    }

    @Test
    public void testReadTwoSampleFiles() {
        Set<File> sampleFiles = new HashSet<>(2);
        sampleFiles.add(new File(getToolTestDataDir(), "samples1.samples"));
        sampleFiles.add(new File(getToolTestDataDir(), "samples2.samples"));
        Collection<String> samples = SampleUtils.getSamplesFromFiles(sampleFiles);
        Assert.assertEquals(samples.size(), 4);
        Assert.assertTrue(samples.contains("NA20845"));
        Assert.assertTrue(samples.contains("NA20846"));
        Assert.assertTrue(samples.contains("NA20847"));
        Assert.assertTrue(samples.contains("NA20849"));
    }

    @Test
    public void testReadOverlappingSampleFiles() {
        Set<File> sampleFiles = new HashSet<>(2);
        sampleFiles.add(new File(getToolTestDataDir(), "samples2.samples"));
        sampleFiles.add(new File(getToolTestDataDir(), "overlapsWithSamples2.samples"));
        Collection<String> samples = SampleUtils.getSamplesFromFiles(sampleFiles);
        Assert.assertEquals(samples.size(), 2);
        Assert.assertTrue(samples.contains("NA20847"));
        Assert.assertTrue(samples.contains("NA20849"));
    }

    @Test
    public void testReadEmptySampleFile() {
        Set<File> sampleFiles = new HashSet<>(1);
        sampleFiles.add(new File(getToolTestDataDir(), "emptySamples.samples"));
        Collection<String> samples = SampleUtils.getSamplesFromFiles(sampleFiles);
        Assert.assertEquals(samples.size(), 0);
    }

    @Test(expectedExceptions=UserException.CouldNotReadInputFile.class)
    public void testReadNonexistentSampleFile() throws Exception {
        Set<File> sampleFiles = new HashSet<>(1);
        sampleFiles.add(BaseTest.getSafeNonExistentFile("nonExistentFile.samples"));
        SampleUtils.getSamplesFromFiles(sampleFiles);
    }
}
