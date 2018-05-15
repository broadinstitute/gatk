package org.broadinstitute.hellbender.tools.spark;

import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Integration tests for {@link PileupSpark}.
 */
public final class PileupSparkIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "walkers/qc/pileup");

    @DataProvider(name = "shuffle")
    public Object[][] shuffleParameters() {
        return new Object[][] { { false }, { true } };
    }

    private File createTempFile() throws IOException {
        final File out = File.createTempFile("out", ".txt");
        out.delete();
        out.deleteOnExit();
        return out;
    }

    @Test(dataProvider = "shuffle")
    public void testSimplePileup(boolean useShuffle) throws Exception {
        final File out = createTempFile();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(NA12878_20_21_WGS_bam);
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(b37_reference_20_21);
        args.add("-L 20:9999900-10000000");
        if (useShuffle) {
            args.add("--shuffle");
            args.add("--num-reducers").add("1");
        }
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_DATA_DIR, "expectedSimplePileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test(dataProvider = "shuffle")
    public void testVerbosePileup(boolean useShuffle) throws Exception {
        final File out = createTempFile();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(NA12878_20_21_WGS_bam);
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(b37_reference_20_21);
        args.add("-L 20:9999990-10000000");
        args.add("-verbose");
        if (useShuffle) {
            args.add("--shuffle");
            args.add("--num-reducers").add("1");
        }
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_DATA_DIR, "expectedVerbosePileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test(dataProvider = "shuffle")
    public void testFeaturesPileup(boolean useShuffle) throws Exception {
        final File out = createTempFile();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(NA12878_20_21_WGS_bam);
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(b37_reference_20_21);
        args.add("-L 20:10000092-10000112");
        args.add("-metadata " + dbsnp_138_b37_20_21_vcf);
        if (useShuffle) {
            args.add("--shuffle");
            args.add("--num-reducers").add("1");
        }
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_DATA_DIR, "expectedFeaturesPileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test(dataProvider = "shuffle")
    public void testInsertLengthPileup(boolean useShuffle) throws Exception {
        final File out = createTempFile();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(NA12878_20_21_WGS_bam);
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(b37_reference_20_21);
        args.add("-L 20:10000092-10000112");
        args.add("--output-insert-length");
        if (useShuffle) {
            args.add("--shuffle");
            args.add("--num-reducers").add("1");
        }
        this.runCommandLine(args.getArgsArray());
        File expected = new File(TEST_DATA_DIR, "expectedInsertLengthPileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test(dataProvider = "shuffle")
    public void testFeaturesPileupHdfs(boolean useShuffle) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path vcfPath = new Path(workingDirectory, "dbsnp_138.b37.20.21.vcf");
            final Path idxPath = new Path(workingDirectory, "dbsnp_138.b37.20.21.vcf.idx");
            cluster.getFileSystem().copyFromLocalFile(new Path(dbsnp_138_b37_20_21_vcf), vcfPath);
            cluster.getFileSystem().copyFromLocalFile(new Path(dbsnp_138_b37_20_21_vcf + ".idx"), idxPath);

            final File out = createTempFile();
            final ArgumentsBuilder args = new ArgumentsBuilder();
            args.add("--input");
            args.add(NA12878_20_21_WGS_bam);
            args.add("--output");
            args.add(out.getAbsolutePath());
            args.add("--reference");
            args.add(b37_reference_20_21);
            args.add("-L 20:10000092-10000112");
            args.add("-metadata " + vcfPath.toString());
            if (useShuffle) {
                args.add("--shuffle");
                args.add("--num-reducers").add("1");
            }
            this.runCommandLine(args.getArgsArray());
            File expected = new File(TEST_DATA_DIR, "expectedFeaturesPileup.txt");
            IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);

        });


    }

}