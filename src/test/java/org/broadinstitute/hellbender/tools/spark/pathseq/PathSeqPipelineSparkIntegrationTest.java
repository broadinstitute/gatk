package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.ValidationStringency;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.annotations.Test;

import java.io.File;

public class PathSeqPipelineSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return PathSeqPipelineSpark.class.getSimpleName();
    }

    public void runTest( final File inputBam, final File expectedBam, final File expectedScores, final boolean isHostAligned) throws Exception {
        final File outputBam = createTempFile("pathseqPipelineTestOutput", ".bam");
        final File outputScores = createTempFile("pathseqPipelineTestOutput", ".txt");

        final String baseResourcePath = "src/test/resources/" + PathSeqBuildKmers.class.getPackage().getName().replace(".", "/");
        final String kmerLibraryPath = baseResourcePath + "/hg19mini.hss";
        final String filterImagePath =  baseResourcePath + "/hg19mini.fasta.bwa_image";

        final File pathogenBwaImage = getTestFile("e_coli_k12_mini.fa.img");
        final File pathogenFasta = getTestFile("e_coli_k12_mini.fa");

        final File taxonomyDatabase = getTestFile("e_coli_k12_mini.db");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(inputBam);
        args.addOutput(outputBam);
        args.addFileArgument("scoresOutputPath", outputScores);
        args.addArgument("kmerLibraryPath", kmerLibraryPath);
        args.addArgument("filterBwaImage", filterImagePath);
        args.addBooleanArgument("isHostAligned", isHostAligned);
        args.addFileArgument("pathogenBwaImage", pathogenBwaImage);
        args.addFileArgument("pathogenFasta", pathogenFasta);
        args.addFileArgument("taxonomicDatabasePath", taxonomyDatabase);

        this.runCommandLine(args);

        SamAssertionUtils.assertEqualBamFiles(outputBam, expectedBam, true, ValidationStringency.STRICT);

        String expectedScoreString = FileUtils.readFileToString(expectedScores);
        String actualScoresString = FileUtils.readFileToString(outputScores);
        PathSeqScoreIntegrationTest.compareScoreTables(expectedScoreString, actualScoresString);
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testPathSeqPipeline() throws Exception {
        final File inputBam = getTestFile("pipeline_input.bam");
        final File expectedBam = getTestFile("pipeline_output.bam");
        final File expectedScores = getTestFile("pipeline_output.txt");
        runTest(inputBam, expectedBam, expectedScores, false);
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testPathSeqPipelineHostAlignedInput() throws Exception {
        final File inputBam = getTestFile("pipeline_input_aligned.bam");
        final File expectedBam = getTestFile("pipeline_output_aligned.bam");
        final File expectedScores = getTestFile("pipeline_output_aligned.txt");
        runTest(inputBam, expectedBam, expectedScores, true);
    }
}