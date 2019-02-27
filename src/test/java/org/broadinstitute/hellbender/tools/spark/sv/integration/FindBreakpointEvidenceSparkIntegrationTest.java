package org.broadinstitute.hellbender.tools.spark.sv.integration;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Integration test on the SV pipeline as it exists right now [2017-03-06]
 */
public class FindBreakpointEvidenceSparkIntegrationTest extends CommandLineProgramTest {

    private static final class FindBreakpointEvidenceSparkIntegrationTestArgs {
        final String expectedAlignedContigsLoc;
        final String bamLoc;
        final String kmerIgnoreListLoc;
        final String alignerRefIndexImgLoc;
        final String outputDir;
        final float bamCoverage;
        final String svEvidenceFilterType;
        final String svGenomeGapsFile;
        final String svGenomeUmapS100File;

        FindBreakpointEvidenceSparkIntegrationTestArgs(final String expectedAlignedContigsLoc,
                                                       final String bamLoc,
                                                       final String kmerIgnoreListLoc,
                                                       final String alignerRefIndexImgLoc,
                                                       final String outputDir,
                                                       final float bamCoverage,
                                                       final String svEvidenceFilterType,
                                                       final String svGenomeGapsFile,
                                                       final String svGenomeUmapS100File) {
            this.expectedAlignedContigsLoc = expectedAlignedContigsLoc;
            this.bamLoc = bamLoc;
            this.kmerIgnoreListLoc = kmerIgnoreListLoc;
            this.alignerRefIndexImgLoc = alignerRefIndexImgLoc;
            this.outputDir = outputDir;
            this.bamCoverage = bamCoverage;
            this.svEvidenceFilterType = svEvidenceFilterType;
            this.svGenomeGapsFile = svGenomeGapsFile;
            this.svGenomeUmapS100File = svGenomeUmapS100File;
        }

        String getCommandLine() {
            return  " -I " + bamLoc +
                    " -O "                    + outputDir + "/assemblies.sam" +
                    " --aligner-index-image " + alignerRefIndexImgLoc +
                    " --kmers-to-ignore " + kmerIgnoreListLoc +
                    " --breakpoint-intervals " + outputDir + "/intervals" +
                    " --fastq-dir "            + outputDir + "/fastq" +
                    " --target-link-file "      + outputDir + "/targetLinks.bedpe" +
                    " --min-evidence-coverage-ratio " + 15 / bamCoverage +
                    " --min-coherent-evidence-coverage-ratio " + 7 / bamCoverage +
                    " --sv-evidence-filter-type " + svEvidenceFilterType +
                    (svGenomeGapsFile.isEmpty() ? "" : " --sv-genome-gaps-file " + svGenomeGapsFile) +
                    (svGenomeUmapS100File.isEmpty() ? "" : " --sv-genome-umap-s100-file " + svGenomeUmapS100File);
        }

        @Override
        public String toString() {
            return "FindBreakpointEvidenceSparkIntegrationTestArgs{" +
                    "bam-loc='" + bamLoc + '\'' +
                    ", kmer-ignore-list-loc='" + kmerIgnoreListLoc + '\'' +
                    ", aligner-ref-index-img-loc='" + alignerRefIndexImgLoc + '\'' +
                    ", output-dir='" + outputDir + '\'' +
                    ", bam-coverage='" + bamCoverage + '\'' +
                    ", sv-evidence-filter-type='" + svEvidenceFilterType + '\'' +
                    '}';
        }
    }

    @DataProvider(name = "findBreakpointEvidenceSparkIntegrationTest")
    public Object[][] createTestData() throws IOException {

        List<Object[]> tests = new ArrayList<>();
        final File tempDirNew = BaseTest.createTempDir("forNew");
        tempDirNew.deleteOnExit();
        Files.createDirectories(Paths.get(tempDirNew.getAbsolutePath()+"/fastq"));
        // Test pipeline with BreakpointDensityFilter
        tests.add(
                new Object[]{
                        new FindBreakpointEvidenceSparkIntegrationTestArgs(
                                SVIntegrationTestDataProvider.TEST_CONTIG_SAM,
                                SVIntegrationTestDataProvider.TEST_BAM, SVIntegrationTestDataProvider.KMER_KILL_LIST,
                                SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG, tempDirNew.getAbsolutePath(),
                                SVIntegrationTestDataProvider.TEST_BAM_COVERAGE, SVIntegrationTestDataProvider.DENSITY_FILTER,
                                "", ""
                        )
                }
        );
        // Test pipeline using classifier filter
        tests.add(
                new Object[]{
                        new FindBreakpointEvidenceSparkIntegrationTestArgs(
                                SVIntegrationTestDataProvider.TEST_CONTIG_SAM_CLASSIFIER,
                                SVIntegrationTestDataProvider.TEST_BAM, SVIntegrationTestDataProvider.KMER_KILL_LIST,
                                SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG, tempDirNew.getAbsolutePath(),
                                SVIntegrationTestDataProvider.TEST_BAM_COVERAGE, SVIntegrationTestDataProvider.CLASSIFIER_FILTER,
                                SVIntegrationTestDataProvider.TEST_GENOME_GAPS_FILE, SVIntegrationTestDataProvider.TEST_GENOME_UMAP100_FILE
                        )
                }
        );

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "findBreakpointEvidenceSparkIntegrationTest", groups = "sv")
    public void testFindBreakpointRunnableLocal(final FindBreakpointEvidenceSparkIntegrationTestArgs params) throws IOException {
        final List<String> args = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
        runCommandLine(args);

        final File actualSAMfile = new File(params.outputDir, "assemblies.sam");
        final File expectedSAMfile = IOUtils.getPath(params.expectedAlignedContigsLoc).toFile();
        IntegrationTestSpec.assertEqualTextFiles(actualSAMfile, expectedSAMfile);
    }

    @Test(dataProvider = "findBreakpointEvidenceSparkIntegrationTest", groups = "sv")
    public void testFindBreakpointRunnableMiniCluster(final FindBreakpointEvidenceSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLine()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            // copy local data to minicluster file system and update args
            changeArgCopyFromLocal(argsToBeModified, "-I", new Path(workingDirectory, "hdfs.bam"), cluster);
            changeArgCopyFromLocal(argsToBeModified, "--kmers-to-ignore", new Path(workingDirectory, "dummy.kill.kmers"), cluster);

            // outputs, prefix with hdfs address
            changeArg(argsToBeModified, "-O",
                    new Path(workingDirectory, "assemblies.sam").toUri().toString());
            changeArg(argsToBeModified, "--breakpoint-intervals",
                    new Path(workingDirectory, "intervals").toUri().toString());
            changeArg(argsToBeModified, "--fastq-dir",
                    new Path(workingDirectory, "fastq").toUri().toString());

            runCommandLine(argsToBeModified);

            final File actualSAMfile = BaseTest.createTempFile("assemblies", ".sam");
            cluster.getFileSystem().copyToLocalFile(new Path(workingDirectory, "assemblies.sam"), new Path(actualSAMfile.toURI()));
            final File expectedSAMfile = IOUtils.getPath(params.expectedAlignedContigsLoc).toFile();
            IntegrationTestSpec.assertEqualTextFiles(actualSAMfile, expectedSAMfile);
        });
    }

    static private void changeArg(final List<String> argsToBeModified, final String arg, final String newVal) {
        final int idx = argsToBeModified.indexOf(arg);
        argsToBeModified.set(idx + 1, newVal);
    }

    static private void changeArgCopyFromLocal(final List<String> argsToBeModified, final String arg, final Path newPath,
                                               final MiniDFSCluster cluster) throws Exception {
        final int idx = argsToBeModified.indexOf(arg);
        final File file = new File(argsToBeModified.get(idx + 1));
        cluster.getFileSystem().copyFromLocalFile(new Path(file.toURI()), newPath);
        argsToBeModified.set(idx + 1, newPath.toString());
    }
}
