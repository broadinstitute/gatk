package org.broadinstitute.hellbender.tools.spark;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.util.Collections;


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
                    "--inputGCSPath " + gcpInputPath +
                            " --outputHDFSDirectory " + tempPath;
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

    @Test(groups = {"spark", "bucket"})
    public void testCopyDirectory() throws Exception {
        MiniDFSCluster cluster = null;
        try {
            final Configuration conf = new Configuration();
            // set the minicluster to have a very low block size so that we can test transferring a file in chunks without actually needing to move a big file
            conf.set("dfs.blocksize", "1048576");
            cluster = MiniClusterUtils.getMiniCluster(conf);

            // copy a directory
            final Path tempPath = MiniClusterUtils.getTempPath(cluster, "test", "dir");

            // directory contains two small files named foo.txt and bar.txt
            final String gcpInputPath = getGCPTestInputPath() + "parallel_copy/";
            String args =
                    "--inputGCSPath " + gcpInputPath +
                            " --outputHDFSDirectory " + tempPath;
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

            Assert.assertEquals(filesFound, 2);


            Assert.assertEquals(Utils.calculateFileMD5(new File(tempDir + "/foo.txt")), "d3b07384d113edec49eaa6238ad5ff00");
            Assert.assertEquals(Utils.calculateFileMD5(new File(tempDir + "/bar.txt")), "c157a79031e1c40f85931829bc5fc552");
        } finally {
            MiniClusterUtils.stopCluster(cluster);
        }
    }

}
