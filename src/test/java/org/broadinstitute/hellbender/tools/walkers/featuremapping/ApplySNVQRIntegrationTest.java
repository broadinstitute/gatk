package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.util.LinkedList;
import java.util.List;

public class ApplySNVQRIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;
    private static final double DOUBLE_EPSILON = 0.001;

    private static String testDir = publicTestDir + FlowTestConstants.APPLY_SNVQR_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testApplySNVQRTest");
        final String filename = "snv_apply_snvqr_output.sam";
        final File expectedFile = new File(testDir + "/" + filename);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + filename);

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", publicTestDir + FlowTestConstants.APPLY_SNVQR_DATA_DIR + "/snv_apply_snvqr_input.bam",
                "-L", "chr1:1-10122",
                "--model", testDir + "/model.json",
                "--conf", testDir + "/config.json",
                "--limit-score", "10",
                "--debug-read-name", "30020185_2-UGAv3-182-1989782468",
                "--verbosity", "INFO",
                "--stats-path-prefix", "/tmp/apply_stats"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the output file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "@PG");
        }
    }

    @Test
    public void testNoModel() throws IOException {

        final File outputDir = createTempDir("testApplySNVQRTest");
        final String filename = "snv_apply_snvqr_output_no_model.sam";
        final File expectedFile = new File(testDir + "/" + filename);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + filename);

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", publicTestDir + FlowTestConstants.APPLY_SNVQR_DATA_DIR + "/snv_apply_snvqr_input.bam",
                "-L", "chr1:1-10122",
                "--conf", testDir + "/config_no_model.json",
                "--limit-score", "10",
                "--debug-read-name", "30020185_2-UGAv3-182-1989782468",
                "--verbosity", "INFO",
                "--stats-path-prefix", "/tmp/apply_stats"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the output file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "@PG");
        }
    }

    @Test
    public void testNoModelNoConf() throws IOException {

        final File outputDir = createTempDir("testApplySNVQRTest");
        final String filename = "snv_apply_snvqr_output_no_model_no_conf.sam";
        final File expectedFile = new File(testDir + "/" + filename);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/" + filename);

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", publicTestDir + FlowTestConstants.APPLY_SNVQR_DATA_DIR + "/snv_apply_snvqr_input.bam",
                "-L", "chr1:1-15000",
                "--limit-score", "10",
                "--debug-read-name", "30020185_2-UGAv3-182-1989782468",
                "--verbosity", "INFO",
                "--stats-path-prefix", "/tmp/apply_stats",
                //"--negative-score-override", "0.0001",
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the output file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "@PG");
        }
    }

    @Test
    public void testXGBoostWorking() throws XGBoostError, IOException {

        // load model
        Booster booster = XGBoost.loadModel(testDir + "/model.json");

        // load test data matrix
        DMatrix testData = new DMatrix(testDir + "/train.buffer");

        // predict
        float[][] predicts = booster.predict(testData);

        // for now, its enough that this did not generate an exception and that there is actual an actual array returned
        Assert.assertNotNull(predicts);
        Assert.assertTrue(predicts.length > 0);
        Assert.assertTrue(predicts[0].length > 0);

        // write to a temp file
        savePredicts(predicts);

        // compare to expected results
        final File expectedFile = new File(testDir + "/y_train_pred.csv");
        BufferedReader reader = new BufferedReader(new FileReader(expectedFile));
        reader.readLine(); // skip heading line
        for ( int row = 0 ; row < predicts.length ; row++ ) {
            final String line = reader.readLine();
            Assert.assertNotNull(line);
            final String[] toks = line.split(",");
            Assert.assertEquals(toks.length, predicts[row].length);
            for ( int col = 0 ; col < toks.length ; col++ ) {
                Assert.assertEquals(predicts[row][col], Double.parseDouble(toks[col]), DOUBLE_EPSILON, "row: " + row + ", col: " + col);
            }
        }
        reader.close();;
    }

    private void savePredicts(float[][] predicts) throws IOException {

        final File outputDir = new File("/tmp"); // createTempDir("testApplySNVQRTest");
        final String outputPath = outputDir + "/test_results.csv";
        System.out.println("saving into: " + outputPath);
        PrintWriter pw = new PrintWriter(outputPath);
        List<String> lineElems = new LinkedList<>();
        for ( int col = 0 ; col < predicts[0].length ; col++ ) {
            lineElems.add(Integer.toString(col));
        }
        pw.println(StringUtils.join(lineElems, ","));
        for ( int row = 0 ; row < predicts.length ; row++ ) {
            lineElems.clear();;
            for ( int col = 0 ; col < predicts[row].length ; col++ ) {
                lineElems.add(Double.toString(predicts[row][col]));
            }
            pw.println(StringUtils.join(lineElems, ","));
        }
        pw.close();
    }
}
