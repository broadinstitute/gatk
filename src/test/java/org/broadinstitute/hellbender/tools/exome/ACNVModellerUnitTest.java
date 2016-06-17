package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Log;
import org.apache.commons.io.FileUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.exome.AllelicCNV.*;

/**
 * Unit tests for {@link ACNVModeller}.  Note that behavior of the tests depends on the content of the input files
 * and that changing them may break the tests.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ACNVModellerUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final String SAMPLE_NAME = "test";
    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR, "coverages-for-acnv-modeller.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR, "snps-for-acnv-modeller.tsv");
    private static final File SEGMENT_FILE = new File(TEST_SUB_DIR, "segments-for-acnv-modeller.seg");  //44 segments (some spurious)
    private static final File SEGMENTS_TRUTH_FILE = new File(TEST_SUB_DIR, "segments-truth-for-acnv-modeller.seg");  //30 true segments

    private static final int NUM_SAMPLES = 100;
    private static final int NUM_BURN_IN = 50;

    private static final double INTERVAL_THRESHOLD = 1.; //merge when parameter modes are within 95% HPD interval
    private static final int MAX_SIMILAR_SEGMENT_MERGE_ITERATIONS = 25;

    private static final List<SimpleInterval> SEGMENTS_TRUTH = SegmentUtils.readIntervalsFromSegmentFile(SEGMENTS_TRUTH_FILE);

    //similar-segment merging without refitting may miss some merges, depending on the input test data; we allow up to 2 misses
    private static final int DELTA_NUM_SEGMENTS_WITHOUT_REFITTING = 2;

    private static final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

    @BeforeTest
    private void doBeforeTest() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
    }

    /**
     * Test of similar-segment merging using only copy-ratio data (simulated coverages and segments).
     * Spurious breakpoints have been introduced into the list of true segments; similar-segment merging should
     * remerge segments broken by these breakpoints and reproduce the original list of true segments.
     */
    @Test
    public void testMergeSimilarSegmentsCopyRatio() throws IOException {
        final boolean doRefit = true;

        final String tempDir = publicTestDir + "similar-segment-copy-ratio-test";
        final File tempDirFile = createTempDir(tempDir);

        //load data (coverages and segments)
        final ReadCountCollection coverage = ReadCountCollectionUtils.parse(COVERAGES_FILE);
        final List<AllelicCount> snpCountsDummy =
                Collections.singletonList(new AllelicCount(new SimpleInterval("1", 1, 1), 0, 1));
        final Genome genome = new Genome(coverage, snpCountsDummy, SAMPLE_NAME);
        final SegmentedGenome segmentedGenome = new SegmentedGenome(SEGMENT_FILE, genome);

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller = new ACNVModeller(segmentedGenome, NUM_SAMPLES, NUM_BURN_IN, 10, 0, ctx);
        //check that model is completely fit at construction
        Assert.assertTrue(modeller.isModelFit());
        //perform iterations of similar-segment merging until all similar segments are merged
        int prevNumSegments;
        for (int numIterations = 1; numIterations <= MAX_SIMILAR_SEGMENT_MERGE_ITERATIONS; numIterations++) {
            prevNumSegments = modeller.getACNVModeledSegments().size();
            modeller.performSimilarSegmentMergingIteration(INTERVAL_THRESHOLD, Double.POSITIVE_INFINITY, doRefit);
            if (modeller.getACNVModeledSegments().size() == prevNumSegments) {
                break;
            }
        }
        //write final segments to file
        final File finalModeledSegmentsFile =
                new File(tempDirFile.getAbsolutePath() + "/test-" + FINAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        modeller.writeACNVModeledSegmentFile(finalModeledSegmentsFile);

        //check equality of segments
        final List<SimpleInterval> segmentsResult = SegmentUtils.readIntervalsFromSegmentFile(finalModeledSegmentsFile);
        Assert.assertEquals(segmentsResult, SEGMENTS_TRUTH);
    }

    /**
     * Test of similar-segment merging using simulated data (coverages, SNP counts, and segments).
     * Spurious breakpoints have been introduced into the list of true segments; similar-segment merging should
     * remerge segments broken by these breakpoints and reproduce the original list of true segments.
     */
    @Test
    public void testMergeSimilarSegmentsWithRefitting() {
        final boolean doRefit = true;

        final String tempDir = publicTestDir + "similar-segment-test";
        final File tempDirFile = createTempDir(tempDir);

        //load data (coverages, SNP counts, and segments)
        final Genome genome = new Genome(COVERAGES_FILE, TUMOR_ALLELIC_COUNTS_FILE, SAMPLE_NAME);
        final SegmentedGenome segmentedGenome = new SegmentedGenome(SEGMENT_FILE, genome);

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller =
                new ACNVModeller(segmentedGenome, NUM_SAMPLES, NUM_BURN_IN, NUM_SAMPLES, NUM_BURN_IN, ctx);
        //check that model is completely fit at construction
        Assert.assertTrue(modeller.isModelFit());
        //perform iterations of similar-segment merging until all similar segments are merged
        int prevNumSegments;
        for (int numIterations = 1; numIterations <= MAX_SIMILAR_SEGMENT_MERGE_ITERATIONS; numIterations++) {
            prevNumSegments = modeller.getACNVModeledSegments().size();
            modeller.performSimilarSegmentMergingIteration(INTERVAL_THRESHOLD, INTERVAL_THRESHOLD, doRefit);
            if (modeller.getACNVModeledSegments().size() == prevNumSegments) {
                break;
            }
        }
        //write final segments to file
        final File finalModeledSegmentsFile =
                new File(tempDirFile.getAbsolutePath() + "/test-" + FINAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        modeller.writeACNVModeledSegmentFile(finalModeledSegmentsFile);

        //check equality of segments
        final List<SimpleInterval> segmentsResult = SegmentUtils.readIntervalsFromSegmentFile(finalModeledSegmentsFile);
        final List<SimpleInterval> segmentsTruth = SegmentUtils.readIntervalsFromSegmentFile(SEGMENTS_TRUTH_FILE);
        Assert.assertEquals(segmentsResult, segmentsTruth);
    }

    /**
     * Same as {@link ACNVModellerUnitTest#testMergeSimilarSegmentsWithRefitting()}, but without intermediate model refitting.
     */
    @Test
    public void testMergeSimilarSegmentsWithoutRefitting() {
        final boolean doRefit = false;

        final String tempDir = publicTestDir + "similar-segment-test-without-refitting";
        final File tempDirFile = createTempDir(tempDir);

        //load data (coverages, SNP counts, and segments)
        final Genome genome = new Genome(COVERAGES_FILE, TUMOR_ALLELIC_COUNTS_FILE, SAMPLE_NAME);
        final SegmentedGenome segmentedGenome = new SegmentedGenome(SEGMENT_FILE, genome);

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller =
                new ACNVModeller(segmentedGenome, NUM_SAMPLES, NUM_BURN_IN, NUM_SAMPLES, NUM_BURN_IN, ctx);
        //check that model is completely fit at construction
        Assert.assertTrue(modeller.isModelFit());
        //perform iterations of similar-segment merging until all similar segments are merged
        int prevNumSegments;
        for (int numIterations = 1; numIterations <= MAX_SIMILAR_SEGMENT_MERGE_ITERATIONS; numIterations++) {
            prevNumSegments = modeller.getACNVModeledSegments().size();
            modeller.performSimilarSegmentMergingIteration(INTERVAL_THRESHOLD, INTERVAL_THRESHOLD, doRefit);
            if (modeller.getACNVModeledSegments().size() == prevNumSegments) {
                break;
            }
        }
        //check that model is not completely fit at this point
        Assert.assertTrue(!modeller.isModelFit());

        //write final segments to file; this should force model refit
        final File finalModeledSegmentsFile =
                new File(tempDirFile.getAbsolutePath() + "/test-" + FINAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        modeller.writeACNVModeledSegmentFile(finalModeledSegmentsFile);
        Assert.assertTrue(modeller.isModelFit());

        //check final number of segments (we may miss some merges by not performing intermediate refits)
        final List<SimpleInterval> segmentsResult = SegmentUtils.readIntervalsFromSegmentFile(finalModeledSegmentsFile);
        Assert.assertEquals(segmentsResult.size(), SEGMENTS_TRUTH.size(), DELTA_NUM_SEGMENTS_WITHOUT_REFITTING);
    }

    /**
     * Test for {@link ACNVModeller#writeACNVModeledSegmentFile} and {@link ACNVModeller#writeACNVModelParameterFiles}.
     */
    @Test
    public void testSegmentAndParameterFileOutput() {
        final String tempDir = publicTestDir + "acnv-modeller-file-output";
        final File tempDirFile = createTempDir(tempDir);

        //load data (coverages, SNP counts, and segments)
        final Genome genome = new Genome(COVERAGES_FILE, TUMOR_ALLELIC_COUNTS_FILE, SAMPLE_NAME);
        final SegmentedGenome segmentedGenome = new SegmentedGenome(SEGMENT_FILE, genome);

        //just do a few iterations
        final int numSamples = 10;
        final int numBurnIn = 0;

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller =
                new ACNVModeller(segmentedGenome, numSamples, numBurnIn, numSamples, numBurnIn, ctx);

        //write segments to file
        final File modeledSegmentsFile = new File(tempDirFile.getAbsolutePath() + "/test-" + FINAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        modeller.writeACNVModeledSegmentFile(modeledSegmentsFile);

        //write parameters to file
        final File copyRatioParameterFile = new File(tempDirFile.getAbsolutePath() + "/test-" + FINAL_FIT_FILE_TAG + CR_PARAMETER_FILE_SUFFIX);
        final File alleleFractionParameterFile = new File(tempDirFile.getAbsolutePath() + "/test-" + FINAL_FIT_FILE_TAG + AF_PARAMETER_FILE_SUFFIX);
        modeller.writeACNVModelParameterFiles(copyRatioParameterFile, alleleFractionParameterFile);

        //check that all files are created with appropriate size
        Assert.assertTrue(modeledSegmentsFile.isFile(), modeledSegmentsFile.getAbsolutePath() + " is not a file.");
        Assert.assertTrue(copyRatioParameterFile.isFile(), copyRatioParameterFile.getAbsolutePath() + " is not a file.");
        Assert.assertTrue(alleleFractionParameterFile.isFile(), alleleFractionParameterFile.getAbsolutePath() + " is not a file.");
        try {
            Assert.assertTrue(FileUtils.readLines(modeledSegmentsFile).size() >= 2);            //column header + at least one segment
            Assert.assertTrue(FileUtils.readLines(copyRatioParameterFile).size() == 3);         //column header + 2 global parameters
            Assert.assertTrue(FileUtils.readLines(alleleFractionParameterFile).size() == 4);    //column header + 3 global parameters
        } catch (final IOException e) {
            Assert.fail("Error reading file size.", e);
        }
    }
}