package org.broadinstitute.hellbender.utils.dataflow;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.util.IOUtil;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.GeneralSecurityException;

public final class BucketUtilsTest extends BaseTest {

    @Test
    public void testIsCloudStorageURL(){
        Assert.assertTrue(BucketUtils.isCloudStorageUrl("gs://a_bucket/bucket"));
        Assert.assertFalse(BucketUtils.isCloudStorageUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isCloudStorageUrl("localFile"));
    }

    @Test
    public void testIsHadoopURL(){
        Assert.assertFalse(BucketUtils.isHadoopUrl("gs://a_bucket/bucket"));
        Assert.assertTrue(BucketUtils.isHadoopUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isHadoopUrl("localFile"));
    }

    @Test
    public void testIsRemoteStorageURL(){
        Assert.assertTrue(BucketUtils.isRemoteStorageUrl("gs://a_bucket/bucket"));
        Assert.assertTrue(BucketUtils.isRemoteStorageUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isRemoteStorageUrl("localFile"));
    }

    @Test
    public void testCopyLocal() throws IOException {
        final String src = publicTestDir+"empty.vcf";
        File dest = createTempFile("copy-empty",".vcf");

        BucketUtils.copyFile(src, null, dest.getPath());
        IOUtil.assertFilesEqual(new File(src), dest);
    }

    @Test
    public void testDeleteLocal() throws IOException, GeneralSecurityException {
        File dest = createTempFile("temp-fortest",".txt");
        try (FileWriter fw = new FileWriter(dest)){
            fw.write("Goodbye, cruel world!");
        }
        BucketUtils.deleteFile(dest.getPath(), null);
        if (dest.exists()) {
            Assert.fail("File '"+dest.getPath()+"' was not properly deleted as it should.");
        }
    }

    @Test(groups={"bucket"})
    public void testCopyAndDeleteGCS() throws IOException, GeneralSecurityException {
        final String src = publicTestDir + "empty.vcf";
        File dest = createTempFile("copy-empty", ".vcf");
        final String intermediate = BucketUtils.randomRemotePath(getDataflowTestStaging(), "test-copy-empty", ".vcf");
        Assert.assertTrue(BucketUtils.isCloudStorageUrl(intermediate), "!BucketUtils.isCloudStorageUrl(intermediate)");
        GCSOptions popts = getAuthenticatedPipelineOptions();
        BucketUtils.copyFile(src, popts, intermediate);
        BucketUtils.copyFile(intermediate, popts, dest.getPath());
        IOUtil.assertFilesEqual(new File(src), dest);
        Assert.assertTrue(BucketUtils.fileExists(intermediate, popts));
        BucketUtils.deleteFile(intermediate, popts);
        Assert.assertFalse(BucketUtils.fileExists(intermediate, popts));
    }

    @Test
    public void testCopyAndDeleteHDFS() throws IOException, GeneralSecurityException {
        final String src = publicTestDir + "empty.vcf";
        File dest = createTempFile("copy-empty", ".vcf");

        MiniDFSCluster cluster = null;
        try {
            cluster = new MiniDFSCluster.Builder(new Configuration()).build();
            String staging = cluster.getFileSystem().getWorkingDirectory().toString();
            final String intermediate = BucketUtils.randomRemotePath(staging, "test-copy-empty", ".vcf");
            Assert.assertTrue(BucketUtils.isHadoopUrl(intermediate), "!BucketUtils.isHadoopUrl(intermediate)");

            PipelineOptions popts = null;
            BucketUtils.copyFile(src, popts, intermediate);
            BucketUtils.copyFile(intermediate, popts, dest.getPath());
            IOUtil.assertFilesEqual(new File(src), dest);
            Assert.assertTrue(BucketUtils.fileExists(intermediate, popts));
            BucketUtils.deleteFile(intermediate, popts);
            Assert.assertFalse(BucketUtils.fileExists(intermediate, popts));
        } finally {
            if (cluster != null) {
                cluster.shutdown();
            }
        }
    }
}