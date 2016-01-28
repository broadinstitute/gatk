package org.broadinstitute.hellbender.utils.gcs;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
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
    public void testIsFileURL(){
        Assert.assertTrue(BucketUtils.isFileUrl("file:///somefile/something"));
        Assert.assertTrue(BucketUtils.isFileUrl("file:/something"));
        Assert.assertFalse(BucketUtils.isFileUrl("gs://a_bucket"));
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
        Assert.assertFalse(dest.exists(), "File '"+dest.getPath()+"' was not properly deleted as it should.");
    }

    @Test(groups={"bucket"})
    public void testCopyAndDeleteGCS() throws IOException, GeneralSecurityException {
        final String src = publicTestDir + "empty.vcf";
        File dest = createTempFile("copy-empty", ".vcf");
        final String intermediate = BucketUtils.randomRemotePath(getGCPTestStaging(), "test-copy-empty", ".vcf");
        Assert.assertTrue(BucketUtils.isCloudStorageUrl(intermediate), "!BucketUtils.isCloudStorageUrl(intermediate)");
        PipelineOptions popts = getAuthenticatedPipelineOptions();
        BucketUtils.copyFile(src, popts, intermediate);
        BucketUtils.copyFile(intermediate, popts, dest.getPath());
        IOUtil.assertFilesEqual(new File(src), dest);
        Assert.assertTrue(BucketUtils.fileExists(intermediate, popts));
        BucketUtils.deleteFile(intermediate, popts);
        Assert.assertFalse(BucketUtils.fileExists(intermediate, popts));
    }

    @Test
    public void testCopyAndDeleteHDFS() throws Exception {
        final String src = publicTestDir + "empty.vcf";
        File dest = createTempFile("copy-empty", ".vcf");

        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final String intermediate = BucketUtils.randomRemotePath(MiniClusterUtils.getWorkingDir(cluster).toString(), "test-copy-empty", ".vcf");
            Assert.assertTrue(BucketUtils.isHadoopUrl(intermediate), "!BucketUtils.isHadoopUrl(intermediate)");

            PipelineOptions popts = null;
            BucketUtils.copyFile(src, popts, intermediate);
            BucketUtils.copyFile(intermediate, popts, dest.getPath());
            IOUtil.assertFilesEqual(new File(src), dest);
            Assert.assertTrue(BucketUtils.fileExists(intermediate, popts));
            BucketUtils.deleteFile(intermediate, popts);
            Assert.assertFalse(BucketUtils.fileExists(intermediate, popts));
        });
    }

    @Test
    public void testDirSize() throws IOException {
        File dir = createTempDir("dir");
        File file1 = new File(dir, "file1.txt");
        File file2 = new File(dir, "file2.txt");
        File subdir = new File(dir, "sub");
        subdir.mkdir();
        File file3 = new File(subdir, "file3.txt");

        for (File file : new File[] { file1, file2, file3 }) {
            try (FileWriter fw = new FileWriter(file)){
                fw.write("Hello!");
            }
        }

        long fileSize = BucketUtils.fileSize(file1.getAbsolutePath(), null);
        Assert.assertTrue(fileSize > 0);
        long dirSize = BucketUtils.dirSize(dir.getAbsolutePath(), null);
        Assert.assertEquals(dirSize, fileSize * 2);
    }

}
