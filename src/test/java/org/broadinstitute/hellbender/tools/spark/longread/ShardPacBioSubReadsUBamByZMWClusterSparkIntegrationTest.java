package org.broadinstitute.hellbender.tools.spark.longread;

import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.DistributedFileSystem;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class ShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTest extends CommandLineProgramTest {

    private static final String THIS_TEST_FOLDER = largeFileTestDir + "longreads/";

    private static final String pbCCSUnalignedReads = THIS_TEST_FOLDER +
            "NA12878rep1.chrM_and_chr20.tiny_subset.unaligned.subreads.bam";

    private static final String pbCCSUnalignedReadsBAI = THIS_TEST_FOLDER +
            "NA12878rep1.chrM_and_chr20.tiny_subset.unaligned.subreads.bam.bai";

    private static final String pbCCSUnalignedReadsSBI = THIS_TEST_FOLDER +
            "NA12878rep1.chrM_and_chr20.tiny_subset.unaligned.subreads.bam.sbi";

    private static final String OUTPUT_SPILT_PREFIX = "PBCCSUBAMSparkSplit";

    private static final class ShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTestArgs {
        final String workDir;

        ShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTestArgs(final String workDir) {
            this.workDir = workDir;
        }

        String getCommandLine() {
            return  " -I " + pbCCSUnalignedReads +
                    " --read-index " + pbCCSUnalignedReadsSBI +
                    " -O " + workDir + "/" + OUTPUT_SPILT_PREFIX;
        }
    }

    @DataProvider(name = "forShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTest")
    private Object[][] createData() {
        List<Object[]> data = new ArrayList<>();

        final File testWorkDir = createTempDir("ShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTest");
        final ShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTestArgs testArgs =
                new ShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTestArgs(testWorkDir.getAbsolutePath());
        data.add(new Object[]{testArgs});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "longreads", dataProvider = "forShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTest")
    public void testRunLocal(final ShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTestArgs params) throws Exception {
        final List<String> args = Arrays.asList( new ArgumentsBuilder().addRaw(params.getCommandLine()).getArgsArray() );
        runCommandLine(args);

        final List<java.nio.file.Path> splits = Files.list(Paths.get(params.workDir)).collect(Collectors.toList());
        testReadsOfSameZMWAreNotSeparated(Paths.get(pbCCSUnalignedReads), splits);
    }

    @Test(groups = "longreads", dataProvider = "forShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTest")
    public void testRunHDFS(final ShardPacBioSubReadsUBamByZMWClusterSparkIntegrationTestArgs params) throws Exception {
        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {
            final DistributedFileSystem fileSystem = cluster.getFileSystem();

            final List<String> argsToBeModified = Arrays.asList( new ArgumentsBuilder().addRaw(params.getCommandLine()).getArgsArray() );
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            int idx = 0;

            // copy inputs to HDFS
            idx = argsToBeModified.indexOf("-I");
            Path path = new Path(workingDirectory, "hdfs.bam");
            File bamFile = new File(argsToBeModified.get(idx+1));
            fileSystem.copyFromLocalFile(new Path(bamFile.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            idx = argsToBeModified.indexOf("--read-index");
            path = new Path(workingDirectory, "hdfs.bam.sbi");
            File sbiFile = new File(argsToBeModified.get(idx+1));
            fileSystem.copyFromLocalFile(new Path(sbiFile.toURI()), path);
            argsToBeModified.set(idx+1, path.toUri().toString());

            // outputs, prefix with hdfs address
            idx = argsToBeModified.indexOf("-O");
            Path outputDir = new Path(workingDirectory, "split");
            String outputPrefix = outputDir.toUri().toString() + "/test";
            argsToBeModified.set(idx+1, outputPrefix);

            runCommandLine(argsToBeModified);

            // copy back the split bams from HDFS and test
            final File hdfsToLocal = GATKBaseTest.createTempDir("hdfsToLocal");
            hdfsToLocal.deleteOnExit();
            for(final FileStatus status : fileSystem.listStatus(outputDir)){
                fileSystem.copyToLocalFile(status.getPath(), new Path(hdfsToLocal.toURI()));
            }
            final List<java.nio.file.Path> splits = Files.list(Paths.get(hdfsToLocal.toString())).collect(Collectors.toList());
            testReadsOfSameZMWAreNotSeparated(Paths.get(pbCCSUnalignedReads), splits);
        });
    }

    /**
     * To test that
     * <ul>
     *     <li>no reads are lost, and</li>
     *     <li>no ZMWs have their associated reads dumped into separate BAMs</li>
     * </ul>
     */
    private static void testReadsOfSameZMWAreNotSeparated(final java.nio.file.Path originalReads,
                                                          final List<java.nio.file.Path> splitBAMs) throws Exception {

        // collect a map from ZMW number to the number of associated reads
        final Map<String, Integer> zmwReadCounter = new HashMap<>();
        try ( final ReadsPathDataSource readsSource = new ReadsPathDataSource(Paths.get(pbCCSUnalignedReads))) {
            for (final GATKRead read : readsSource ) {
                final String zmw = read.getAttributeAsString("zm");
                zmwReadCounter.merge(zmw, 1, Integer::sum);
            }
        }

        // iterate through split uBAMs and check on integrity
        for (final java.nio.file.Path path: splitBAMs) {
            if (!path.endsWith("bam")) continue;
            final Map<String, Integer> counterForThisSplit = new HashMap<>();
            try ( final ReadsPathDataSource readsSource = new ReadsPathDataSource(path)) {
                for (final GATKRead read : readsSource ) {
                    final String zmw = read.getAttributeAsString("zm");
                    counterForThisSplit.merge(zmw, 1, Integer::sum);
                }
            }
            for (String zmw : counterForThisSplit.keySet()) {
                Integer expected = zmwReadCounter.get(zmw);
                Integer actual = counterForThisSplit.get(zmw);
                Assert.assertEquals(actual, expected,
                                    String.format("Read count for ZMW %s doesn't match!", zmw));
            }
        }
    }
}
