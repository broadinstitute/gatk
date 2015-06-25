package org.broadinstitute.hellbender.utils.dataflow;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
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
    public void testCopyLocal() throws IOException {
        final String src = publicTestDir+"empty.vcf";
        File dest = createTempFile("copy-empty",".vsf");

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
        File dest = createTempFile("copy-empty", ".vsf");
        final String intermediate = BucketUtils.randomGcsPath(getDataflowTestStaging(), "test-copy-empty", ".vsf");
        Assert.assertTrue(BucketUtils.isCloudStorageUrl(intermediate), "!BucketUtils.isCloudStorageUrl(intermediate)");
        GCSOptions popts = PipelineOptionsFactory.as(GCSOptions.class);
        popts.setApiKey(getDataflowTestApiKey());
        BucketUtils.copyFile(src, popts, intermediate);
        BucketUtils.copyFile(intermediate, popts, dest.getPath());
        IOUtil.assertFilesEqual(new File(src), dest);
        Assert.assertTrue(BucketUtils.fileExists(intermediate, popts));
        BucketUtils.deleteFile(intermediate, popts);
        Assert.assertFalse(BucketUtils.fileExists(intermediate, popts));
    }
}