package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.charset.StandardCharsets;

public class PathSeqPipelineSparkIntegrationTest extends CommandLineProgramTest {

    static final String baseResourcePath = publicTestDir + PathSeqPipelineSpark.class.getPackage().getName().replace(".", "/");
    static final String kmerLibraryPath = baseResourcePath + "/hg19mini.hss";
    static final String filterImagePath =  publicTestDir + "/hg19mini.fasta.img";

    @DataProvider(name = "pathseqPipelineTestData")
    public Object[][] getTestData() {
        return new Object[][]{
                {"pipeline_input.bam",
                 "pipeline_output.bam",
                 "pipeline_output.txt",
                 "pipeline_output.filter.metrics",
                 "pipeline_output.score.metrics",
                 false},
                {"pipeline_input_aligned.bam",
                 "pipeline_output_aligned.bam",
                 "pipeline_output_aligned.txt",
                 "pipeline_output_aligned.filter.metrics",
                 "pipeline_output_aligned.score.metrics",
                 true}
        };
    }

    @Override
    public String getTestedClassName() {
        return PathSeqPipelineSpark.class.getSimpleName();
    }

    @Test(dataProvider = "pathseqPipelineTestData")
    public void testPipelineTool( final String inputBamFilename, final String expectedBamFilename, final String expectedScoresFilename,
                         final String expectedFilterMetricsFilename, final String expectedScoreMetricsFilename,
                          final boolean isHostAligned) throws Exception {

        final File inputBamFile = getTestFile(inputBamFilename);
        final File expectedBamFile = getTestFile(expectedBamFilename);
        final File expectedScoresFile = getTestFile(expectedScoresFilename);
        final File expectedFilterMetricsFile = getTestFile(expectedFilterMetricsFilename);
        final File expectedScoreMetricsFile = getTestFile(expectedScoreMetricsFilename);

        final File outputBamFile = createTempFile("pathseqPipelineTestOutput", ".bam");
        final File outputScoresFile = createTempFile("pathseqPipelineTestOutput", ".txt");
        final File outputFilterMetricsFile = createTempFile("filter", ".metrics");
        final File outputScoreMetricsFile = createTempFile("score", ".metrics");
        final File pathogenBwaImage = getTestFile("e_coli_k12_mini.fa.img");
        final File pathogenFasta = getTestFile("e_coli_k12_mini.fa");
        final File taxonomyDatabase = getTestFile("e_coli_k12_mini.db");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(inputBamFile);
        args.addOutput(outputBamFile);
        args.addFileArgument(PSScoreArgumentCollection.SCORES_OUTPUT_LONG_NAME, outputScoresFile);
        args.addArgument(PSFilterArgumentCollection.KMER_FILE_PATH_LONG_NAME, kmerLibraryPath);
        args.addArgument(PSFilterArgumentCollection.FILTER_BWA_IMAGE_LONG_NAME, filterImagePath);
        args.addBooleanArgument(PSFilterArgumentCollection.IS_HOST_ALIGNED_LONG_NAME, isHostAligned);
        args.addFileArgument(PSBwaArgumentCollection.MICROBE_BWA_IMAGE_LONG_NAME, pathogenBwaImage);
        args.addFileArgument(PSBwaArgumentCollection.MICROBE_FASTA_LONG_NAME, pathogenFasta);
        args.addFileArgument(PSScoreArgumentCollection.TAXONOMIC_DATABASE_LONG_NAME, taxonomyDatabase);
        args.addFileArgument(PSFilterArgumentCollection.FILTER_METRICS_FILE_LONG_NAME, outputFilterMetricsFile);
        args.addFileArgument(PSScoreArgumentCollection.SCORE_METRICS_FILE_LONG_NAME, outputScoreMetricsFile);
        this.runCommandLine(args);

        SamAssertionUtils.assertEqualBamFiles(outputBamFile, expectedBamFile, true, ValidationStringency.STRICT);

        String expectedScoreString = FileUtils.readFileToString(expectedScoresFile, StandardCharsets.UTF_8);
        String actualScoresString = FileUtils.readFileToString(outputScoresFile, StandardCharsets.UTF_8);
        PathSeqScoreIntegrationTest.compareScoreTables(expectedScoreString, actualScoresString);

        Assert.assertTrue(MetricsFile.areMetricsEqual(outputFilterMetricsFile, expectedFilterMetricsFile));
        Assert.assertTrue(MetricsFile.areMetricsEqual(outputScoreMetricsFile, expectedScoreMetricsFile));
    }

}