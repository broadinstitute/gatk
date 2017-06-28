package org.broadinstitute.hellbender.tools.spark.sv.integration;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.sv.sga.ContigCollectionTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AlignAssembledContigsSparkIntegrationTest extends CommandLineProgramTest {

    private static final class AlignAssembledContigsSparkIntegrationTestArgs {
        final String inputAssemblies;
        final String alignmentOutput;
        final String alignerRefIndexImgLoc;

        AlignAssembledContigsSparkIntegrationTestArgs(final String inputAssemblies, final String alignmentOutput, final String alignerRefIndexImgLoc){
            this.inputAssemblies = inputAssemblies;
            this.alignerRefIndexImgLoc = alignerRefIndexImgLoc;
            this.alignmentOutput = alignmentOutput;
        }

        String getCommandLineNoApiKey() {
            return " --inputAssemblyDir " + inputAssemblies +
                    " -O " + alignmentOutput +
                    " --bwamemIndexImage " + alignerRefIndexImgLoc;
        }
    }

    @DataProvider(name = "alignAssembledContigsSparkIntegrationTest")
    public Object[][] createTestData() throws IOException {
        List<Object[]> tests = new ArrayList<>();

        final File tempWorkingDir = BaseTest.createTempDir("alignAssembledContigsSparkIntegrationTest");
        tempWorkingDir.deleteOnExit();

        final File assemblyWithPackedFastaWithLength = Files.createDirectory(Paths.get(tempWorkingDir.getAbsolutePath()+"/"+"assemblyWithPackedFastaWithLength")).toFile();
        try(  final PrintWriter out = new PrintWriter( new File(assemblyWithPackedFastaWithLength, "contents.txt") )  ){
            out.println( ContigCollectionTest.contigsCollectionPackedFastaAndWithAsmId[0][1] );
        }
        tests.add(new Object[]{new AlignAssembledContigsSparkIntegrationTest.AlignAssembledContigsSparkIntegrationTestArgs(assemblyWithPackedFastaWithLength.getAbsolutePath(),
                tempWorkingDir.getAbsolutePath()+"/"+"outputWithLength",
                SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG)});

        final File assemblyWithPackedFastaWithoutLength = Files.createDirectory(Paths.get(tempWorkingDir.getAbsolutePath()+"/"+"assemblyWithPackedFastaWithoutLength")).toFile();
        try(  final PrintWriter out = new PrintWriter( new File(assemblyWithPackedFastaWithoutLength, "contents.txt") )  ){
            out.println( ContigCollectionTest.contigsCollectionPackedFastaAndWithAsmId[0][1] );
        }
        tests.add(new Object[]{new AlignAssembledContigsSparkIntegrationTest.AlignAssembledContigsSparkIntegrationTestArgs(assemblyWithPackedFastaWithoutLength.getAbsolutePath(),
                tempWorkingDir.getAbsolutePath()+"/"+"outputWithoutLength",
                SVIntegrationTestDataProvider.ALIGNER_INDEX_IMG)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "alignAssembledContigsSparkIntegrationTest", groups = "sv")
    public void testAlignAssembledContigsSparkRunnableLocal(final AlignAssembledContigsSparkIntegrationTest.AlignAssembledContigsSparkIntegrationTestArgs params) throws IOException {
        new IntegrationTestSpec(
                new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getString(),
                SVIntegrationTestDataProvider.dummyExpectedFileNames)
                .executeTest("testAlignAssembledContigsSparkRunnableLocal-", this);
    }

    @Test(dataProvider = "alignAssembledContigsSparkIntegrationTest", groups = "sv")
    public void testAlignAssembledContigsSparkRunnableMiniCluster(final AlignAssembledContigsSparkIntegrationTest.AlignAssembledContigsSparkIntegrationTestArgs params) throws Exception {

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().add(params.getCommandLineNoApiKey()).getArgsArray() );
            final org.apache.hadoop.fs.Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;

            // inputs, copy files
            idx = argsToBeModified.indexOf("--inputAssemblyDir");
            org.apache.hadoop.fs.Path path = new org.apache.hadoop.fs.Path(workingDirectory, "contents.txt");
            File file = new File(argsToBeModified.get(idx+1));
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(file.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            path = new org.apache.hadoop.fs.Path(workingDirectory, "reference.fasta");
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(SVIntegrationTestDataProvider.reference.toURI()), path);
            path = new org.apache.hadoop.fs.Path(workingDirectory, "reference.fasta.fai");
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(SVIntegrationTestDataProvider.reference_fai.toURI()), path);
            path = new org.apache.hadoop.fs.Path(workingDirectory, "reference.2bit");
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(SVIntegrationTestDataProvider.reference_2bit.toURI()), path);
            path = new org.apache.hadoop.fs.Path(workingDirectory, "reference.dict");
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(SVIntegrationTestDataProvider.reference_dict.toURI()), path);

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            path = new org.apache.hadoop.fs.Path(workingDirectory, "alignedContigs.txt");
            argsToBeModified.set(idx+1, path.toUri().toString());

            new IntegrationTestSpec(String.join(" ", argsToBeModified), SVIntegrationTestDataProvider.dummyExpectedFileNames)
                    .executeTest("testAlignAssembledContigsSparkRunnableMiniCluster-", this);
        });
    }
}
