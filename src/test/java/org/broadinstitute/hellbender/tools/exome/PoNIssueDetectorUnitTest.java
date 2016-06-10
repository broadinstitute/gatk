package org.broadinstitute.hellbender.tools.exome;

import com.google.common.collect.Sets;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class PoNIssueDetectorUnitTest extends BaseTest {

    private static String TEST_DIR = "src/test/resources/org/broadinstitute/hellbender/tools/exome/";
    // These files were created by adding simulated events to real data.
    private static File TEST_FILE_AMP = new File(TEST_DIR, "events_tn.txt");
    private static File TEST_FILE_DEL = new File(TEST_DIR, "del_events_tn.txt");

    // This file was created manually from real normal sample tn file.
    private static File TEST_NO_SUSPICIOUS_SAMPLES_FILE = new File(TEST_DIR, "no_events_tn_an.txt");

    @Test
    public void testIdentifySamplesWithSuspiciousContigsAmps() {

        final Set<String> gtBlacklistSamples = new HashSet<>();
        gtBlacklistSamples.add("sample_1");
        gtBlacklistSamples.add("sample_2");
        gtBlacklistSamples.add("sample_3");

        ReadCountCollection allCoverageProfiles = null;
        try {
            allCoverageProfiles = ReadCountCollectionUtils.parse(TEST_FILE_AMP);
        } catch (final IOException ioe) {
            Assert.fail("Could not load test file: " + TEST_FILE_AMP, ioe);
        }
        final List<ReadCountCollection> singleSampleTangentNormalizedReadCounts = PoNIssueDetector.createIndividualReadCountCollections(allCoverageProfiles);

        // By the time we are here, input is assumed to have been tangent normalized.
        final List<String> blacklistSamples = PoNIssueDetector.identifySamplesWithSuspiciousContigs(singleSampleTangentNormalizedReadCounts, PoNIssueDetector.getContigToMedianCRMap(allCoverageProfiles));

        final Set<String> resultSamples = new HashSet<>(blacklistSamples);

        Assert.assertEquals(resultSamples.size(), gtBlacklistSamples.size());
        Assert.assertEquals(Sets.difference(resultSamples, gtBlacklistSamples).size(), 0);
    }

    @Test
    public void testIdentifySamplesWithSuspiciousContigsDels() {

        final Set<String> gtBlacklistSamples = new HashSet<>();
        gtBlacklistSamples.add("sample_1");
        gtBlacklistSamples.add("sample_2");
        gtBlacklistSamples.add("sample_3");

        ReadCountCollection allCoverageProfiles = null;
        try {
            allCoverageProfiles = ReadCountCollectionUtils.parse(TEST_FILE_DEL);
        } catch (final IOException ioe) {
            Assert.fail("Could not load test file: " + TEST_FILE_DEL, ioe);
        }
        final List<ReadCountCollection> singleSampleTangentNormalizedReadCounts = PoNIssueDetector.createIndividualReadCountCollections(allCoverageProfiles);

        // By the time we are here, input is assumed to have been tangent normalized.
        final List<String> blacklistSamples = PoNIssueDetector.identifySamplesWithSuspiciousContigs(singleSampleTangentNormalizedReadCounts, PoNIssueDetector.getContigToMedianCRMap(allCoverageProfiles));

        final Set<String> resultSamples = new HashSet<>(blacklistSamples);

        Assert.assertEquals(resultSamples.size(), gtBlacklistSamples.size());
        Assert.assertEquals(Sets.difference(resultSamples, gtBlacklistSamples).size(), 0);
    }

    @Test
    public void testIdentifySamplesWithSuspiciousContigsNoSuspiciousSamples() {
        ReadCountCollection allCoverageProfiles = null;
        try {
            allCoverageProfiles = ReadCountCollectionUtils.parse(TEST_NO_SUSPICIOUS_SAMPLES_FILE);
        } catch (final IOException ioe) {
            Assert.fail("Could not load test file: " + TEST_NO_SUSPICIOUS_SAMPLES_FILE, ioe);
        }
        final List<ReadCountCollection> singleSampleTangentNormalizedReadCounts = PoNIssueDetector.createIndividualReadCountCollections(allCoverageProfiles);

        // By the time we are here, input is assumed to have been tangent normalized.
        final List<String> blacklistSamples = PoNIssueDetector.identifySamplesWithSuspiciousContigs(singleSampleTangentNormalizedReadCounts, PoNIssueDetector.getContigToMedianCRMap(allCoverageProfiles));
        Assert.assertEquals(blacklistSamples.size(), 0);
    }

    @Test
    public void testIdentifySamplesWithSuspiciousContigsAmpsWithSpark() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final Set<String> gtBlacklistSamples = new HashSet<>();
        gtBlacklistSamples.add("sample_1");
        gtBlacklistSamples.add("sample_2");
        gtBlacklistSamples.add("sample_3");

        ReadCountCollection allCoverageProfiles = null;
        try {
            allCoverageProfiles = ReadCountCollectionUtils.parse(TEST_FILE_AMP);
        } catch (final IOException ioe) {
            Assert.fail("Could not load test file: " + TEST_FILE_AMP, ioe);
        }
        final JavaRDD<ReadCountCollection> allSampleTangentNormalizedReadCounts = PoNIssueDetector.createParallelIndividualReadCountCollections(allCoverageProfiles, ctx);

        // By the time we are here, input is assumed to have been tangent normalized.
        final List<String> blacklistSamples = PoNIssueDetector.identifySamplesWithSuspiciousContigs(allSampleTangentNormalizedReadCounts, ctx, PoNIssueDetector.getContigToMedianCRMap(allCoverageProfiles));

        final Set<String> resultSamples = new HashSet<>(blacklistSamples);

        Assert.assertEquals(resultSamples.size(), gtBlacklistSamples.size());
        Assert.assertEquals(Sets.difference(resultSamples, gtBlacklistSamples).size(), 0);
    }

    @Test
    public void testIdentifySamplesWithSuspiciousContigsDelsWithSpark() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final Set<String> gtBlacklistSamples = new HashSet<>();
        gtBlacklistSamples.add("sample_1");
        gtBlacklistSamples.add("sample_2");
        gtBlacklistSamples.add("sample_3");

        ReadCountCollection allCoverageProfiles = null;
        try {
            allCoverageProfiles = ReadCountCollectionUtils.parse(TEST_FILE_DEL);
        } catch (final IOException ioe) {
            Assert.fail("Could not load test file: " + TEST_FILE_DEL, ioe);
        }
        final JavaRDD<ReadCountCollection> allSampleTangentNormalizedReadCounts = PoNIssueDetector.createParallelIndividualReadCountCollections(allCoverageProfiles, ctx);

        // By the time we are here, input is assumed to have been tangent normalized.
        final List<String> blacklistSamples = PoNIssueDetector.identifySamplesWithSuspiciousContigs(allSampleTangentNormalizedReadCounts, ctx, PoNIssueDetector.getContigToMedianCRMap(allCoverageProfiles));

        final Set<String> resultSamples = new HashSet<>(blacklistSamples);

        Assert.assertEquals(resultSamples.size(), gtBlacklistSamples.size());
        Assert.assertEquals(Sets.difference(resultSamples, gtBlacklistSamples).size(), 0);
    }

    @Test
    public void testIdentifySamplesWithSuspiciousContigsNoSuspiciousSamplesWithSpark() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadCountCollection allCoverageProfiles = null;
        try {
            allCoverageProfiles = ReadCountCollectionUtils.parse(TEST_NO_SUSPICIOUS_SAMPLES_FILE);
        } catch (final IOException ioe) {
            Assert.fail("Could not load test file: " + TEST_NO_SUSPICIOUS_SAMPLES_FILE, ioe);
        }
        final JavaRDD<ReadCountCollection> allSampleTangentNormalizedReadCounts = PoNIssueDetector.createParallelIndividualReadCountCollections(allCoverageProfiles, ctx);

        // By the time we are here, input is assumed to have been tangent normalized.
        final List<String> blacklistSamples = PoNIssueDetector.identifySamplesWithSuspiciousContigs(allSampleTangentNormalizedReadCounts, ctx, PoNIssueDetector.getContigToMedianCRMap(allCoverageProfiles));
        Assert.assertEquals(blacklistSamples.size(), 0);
    }
}
