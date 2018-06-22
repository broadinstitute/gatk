package org.broadinstitute.hellbender.tools.spark;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


public class ParallelCopyGCSDirectoryIntoHDFSSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedToolName() {
        return ParallelCopyGCSDirectoryIntoHDFSSpark.class.getSimpleName();
    }

    @Test(groups = {"spark", "bucket"})
    public void testCopyLargeFile() throws Exception {
        MiniDFSCluster cluster = null;
        try {
            final Configuration conf = new Configuration();
            // set the minicluster to have a very low block size so that we can test transferring a file in chunks without actually needing to move a big file
            conf.set("dfs.blocksize", "1048576");
            cluster = MiniClusterUtils.getMiniCluster(conf);

            // copy a multi-block file
            final Path tempPath = MiniClusterUtils.getTempPath(cluster, "test", "dir");
            final String gcpInputPath = getGCPTestInputPath() + "huge/CEUTrio.HiSeq.WGS.b37.NA12878.chr1_4.bam.bai";
            String args =
                    "--" + ParallelCopyGCSDirectoryIntoHDFSSpark.INPUT_GCS_PATH_LONG_NAME + " " + gcpInputPath +
                            " --" + ParallelCopyGCSDirectoryIntoHDFSSpark.OUTPUT_HDFS_DIRECTORY_LONG_NAME + " " + tempPath;
            ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
            IntegrationTestSpec spec = new IntegrationTestSpec(
                    ab.getString(),
                    Collections.emptyList());
            spec.executeTest("testCopyLargeFile-" + args, this);

            final long fileSizeOnGCS = Files.size(IOUtils.getPath(gcpInputPath));


            final String hdfsPath = tempPath + "/" + "CEUTrio.HiSeq.WGS.b37.NA12878.chr1_4.bam.bai";

            org.apache.hadoop.fs.Path outputHdfsDirectoryPath = new org.apache.hadoop.fs.Path(tempPath.toUri());

            try(FileSystem fs = outputHdfsDirectoryPath.getFileSystem(conf)) {
                long chunkSize = ParallelCopyGCSDirectoryIntoHDFSSpark.getChunkSize(fs);
                Assert.assertTrue(fileSizeOnGCS > chunkSize);
            }

            Assert.assertEquals(BucketUtils.fileSize(hdfsPath),
                    fileSizeOnGCS);

            final File tempDir = createTempDir("ParallelCopy");

            BucketUtils.copyFile(hdfsPath, tempDir + "fileFromHDFS.bam.bai");
            Assert.assertEquals(Utils.calculateFileMD5(new File(tempDir + "fileFromHDFS.bam.bai")), "1a6baa5332e98ef1358ac0fb36f46aaf");
        } finally {
            MiniClusterUtils.stopCluster(cluster);
        }
    }

    @DataProvider(name = "directoryCopy")
    public Object[][] getDirectoryParams() {
        final String gcpInputPath = getGCPTestInputPath() + "parallel_copy/";
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{gcpInputPath, null, new String[] { "foo.txt", "bar.txt"}, new String[] { "d3b07384d113edec49eaa6238ad5ff00", "c157a79031e1c40f85931829bc5fc552"}});
        tests.add(new Object[]{gcpInputPath, "foo*", new String[] { "foo.txt" }, new String[] { "d3b07384d113edec49eaa6238ad5ff00" }});
        return tests.toArray(new Object[][]{});
    }


    @Test(groups = {"spark", "bucket"}, dataProvider = "directoryCopy")
    public void testCopyDirectory(final String gcpInputPath,
                                  final String glob,
                                  final String[] expectedFilesCopied,
                                  final String[] expectedMD5s) throws Exception {
        MiniDFSCluster cluster = null;
        try {
            final Configuration conf = new Configuration();
            // set the minicluster to have a very low block size so that we can test transferring a file in chunks without actually needing to move a big file
            conf.set("dfs.blocksize", "1048576");
            cluster = MiniClusterUtils.getMiniCluster(conf);

            // copy a directory
            final Path tempPath = MiniClusterUtils.getTempPath(cluster, "test", "dir");

            // directory contains two small files named foo.txt and bar.txt

            String args =
                    "--" + ParallelCopyGCSDirectoryIntoHDFSSpark.INPUT_GCS_PATH_LONG_NAME + " " + gcpInputPath +
                            (glob == null ? "" : " --" + ParallelCopyGCSDirectoryIntoHDFSSpark.INPUT_GLOB + " " + glob) +
                            " --" + ParallelCopyGCSDirectoryIntoHDFSSpark.OUTPUT_HDFS_DIRECTORY_LONG_NAME + " " + tempPath;
            ArgumentsBuilder ab = new ArgumentsBuilder().add(args);
            IntegrationTestSpec spec = new IntegrationTestSpec(
                    ab.getString(),
                    Collections.emptyList());
            spec.executeTest("testCopyDirectory-" + args, this);

            org.apache.hadoop.fs.Path outputHdfsDirectoryPath = new org.apache.hadoop.fs.Path(tempPath.toUri());

            final File tempDir = createTempDir("ParallelCopyDir");

            int filesFound = 0;
            try(FileSystem fs = outputHdfsDirectoryPath.getFileSystem(conf)) {
                final RemoteIterator<LocatedFileStatus> hdfsCopies = fs.listFiles(outputHdfsDirectoryPath, false);
                while (hdfsCopies.hasNext()) {
                    final FileStatus next =  hdfsCopies.next();
                    final Path path = next.getPath();
                    BucketUtils.copyFile(path.toString(), tempDir + "/" + path.getName());
                    filesFound ++;
                }
            }

            Assert.assertEquals(filesFound, expectedFilesCopied.length);

            for (int i = 0; i < expectedFilesCopied.length; i++) {
                String fileName = expectedFilesCopied[i];
                String md5 = expectedMD5s[i];
                Assert.assertEquals(Utils.calculateFileMD5(new File(tempDir + "/" + fileName)), md5);
            }



        } finally {
            MiniClusterUtils.stopCluster(cluster);
        }
    }

}
