package org.broadinstitute.hellbender.utils.gcs;

import com.google.cloud.storage.StorageException;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider;
import com.google.cloud.storage.contrib.nio.SeekableByteChannelPrefetcher;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.GeneralSecurityException;
import java.util.function.Function;
import java.util.stream.Stream;

public final class BucketUtilsUnitTest extends GATKBaseTest {

    public static final String LINE_IN_REMOTE_TEXT_FILE = "The Project Gutenberg EBook of The Adventures of Sherlock Holmes";

    /**
     * This file is in a public requester pays bucket and is owned by the broad-gatk-test project.  It must be owned by
     * a different project than the service account doing the testing or the test may fail because it can access the
     * file directly    through alternative permissions.
     */
    public static final String FILE_IN_REQUESTER_PAYS_BUCKET = getGCPRequesterPaysBucket() + "test/resources/nio/big.txt";

    static {
        setDefaultNioOptions();
    }

    private static void setDefaultNioOptions() {
        BucketUtils.setGlobalNIODefaultOptions(
            ConfigFactory.getInstance().getGATKConfig().gcsMaxRetries(),
            ConfigFactory.getInstance().getGATKConfig().gcsProjectForRequesterPays());
    }

    private static void setNoProjectForRequesterPays() {
        BucketUtils.setGlobalNIODefaultOptions(
                ConfigFactory.getInstance().getGATKConfig().gcsMaxRetries(),
                null);
    }

    private static void setRequesterPays(){
        BucketUtils.setGlobalNIODefaultOptions(
                ConfigFactory.getInstance().getGATKConfig().gcsMaxRetries(),
                BaseTest.getGCPTestProject());
    }

    @Test(groups={"bucket"})
    public void testIsGcsUrl(){
        Assert.assertTrue(BucketUtils.isGcsUrl("gs://abucket/bucket"));
        Assert.assertFalse(BucketUtils.isGcsUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isGcsUrl("localFile"));
        Assert.assertFalse(BucketUtils.isGcsUrl("https://www.somewhere.com"));
        Assert.assertFalse(BucketUtils.isGcsUrl("http://www.somewhere.com"));
    }

    @Test
    public void testIsHttpUrl(){
        Assert.assertFalse(BucketUtils.isHttpUrl("gs://abucket/bucket"));
        Assert.assertFalse(BucketUtils.isHttpUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isHttpUrl("localFile"));
        Assert.assertTrue(BucketUtils.isHttpUrl("https://www.somewhere.com"));
        Assert.assertTrue(BucketUtils.isHttpUrl("http://www.somewhere.com"));
    }
    @DataProvider
    public Object[][] getVariousPathsForPrefetching(){
        return new Object[][]{
                {"localFile", false},
                {"file:///local/file", false},
                {"http://www.somewhere.com", true},
                {"https://www.somewhere.com", true},
                {"gs://abucket/bucket", true}
        };
    }

    @Test(groups="bucket", dataProvider = "getVariousPathsForPrefetching")
    public void testIsEligibleForPrefetching(String path, boolean isPrefetchable){
        Assert.assertEquals(BucketUtils.isEligibleForPrefetching(IOUtils.getPath(path)), isPrefetchable);
    }

    @Test(groups="bucket", dataProvider = "getVariousPathsForPrefetching")
    public void testIsEligibleForPrefetchingWithPathSpecifier(String path, boolean isPrefetchable){
        Assert.assertEquals(BucketUtils.isEligibleForPrefetching(new GATKPath(path)), isPrefetchable);
    }

    @Test
    public void testIsEligibleForPrefetchingHdfs(){
        Assert.assertFalse(BucketUtils.isEligibleForPrefetching(new GATKPath("hdfs://namenode/path/to/file")));
    }

    @Test
    public void testGetCloudStorageConfiguration() {
        String mockProject = "yes";
        int mockReopens = 100;
        CloudStorageConfiguration config = BucketUtils.getCloudStorageConfiguration(mockReopens, mockProject);
        Assert.assertEquals(config.maxChannelReopens(), mockReopens);
        Assert.assertEquals(config.userProject(), mockProject);
    }

    @DataProvider
    public Object[][] getPathsForAbsoluteness(){
        return new Object[][]{
                {"some/relative/path", false},
                {"file:///someAbsolutePath", true},
                {"https://some.relative.path.com", true}
        };
    }

    @Test(dataProvider = "getPathsForAbsoluteness")
    public void testMakeFilePathAbsolute(String path, boolean isAbsolute){
        final String absolutized = BucketUtils.makeFilePathAbsolute(path);
        Assert.assertEquals(absolutized.equals(path), isAbsolute, absolutized + " vs " + path);
    }

    @Test
    public void testIsHadoopURL(){
        Assert.assertFalse(BucketUtils.isHadoopUrl("gs://abucket/bucket"));
        Assert.assertTrue(BucketUtils.isHadoopUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isHadoopUrl("localFile"));
        Assert.assertFalse(BucketUtils.isHadoopUrl("http://www.somewhere.com"));
        Assert.assertFalse(BucketUtils.isHadoopUrl("https://www.somewhere.com"));
    }

    @Test
    public void testIsRemoteStorageURL(){
        Assert.assertTrue(BucketUtils.isRemoteStorageUrl("gs://abucket/bucket"));
        Assert.assertTrue(BucketUtils.isRemoteStorageUrl("hdfs://namenode/path/to/file"));
        Assert.assertFalse(BucketUtils.isRemoteStorageUrl("localFile"));
        Assert.assertFalse(BucketUtils.isRemoteStorageUrl("file:///absolutelylocal/file"));
        Assert.assertTrue(BucketUtils.isRemoteStorageUrl("http://www.somewhere.com"));
        Assert.assertTrue(BucketUtils.isRemoteStorageUrl("https://www.somewhere.com"));
    }

    @Test
    public void testIsFileURL(){
        Assert.assertTrue(BucketUtils.isFileUrl("file:///somefile/something"));
        Assert.assertTrue(BucketUtils.isFileUrl("file:/something"));
        Assert.assertFalse(BucketUtils.isFileUrl("gs://abucket"));
        Assert.assertFalse(BucketUtils.isFileUrl("http://www.somewhere.com"));
        Assert.assertFalse(BucketUtils.isFileUrl("https://www.somewhere.com"));
    }

    @Test
    public void testCopyLocal() throws IOException {
        final String src = publicTestDir+"empty.vcf";
        File dest = createTempFile("copy-empty",".vcf");

        BucketUtils.copyFile(src, dest.getPath());
        IOUtil.assertFilesEqual(new File(src), dest);
    }

    @Test
    public void testDeleteLocal() throws IOException {
        File dest = createTempFile("temp-fortest",".txt");
        try (FileWriter fw = new FileWriter(dest)){
            fw.write("Goodbye, cruel world!");
        }
        BucketUtils.deleteFile(dest.getPath());
        Assert.assertFalse(dest.exists(), "File '"+dest.getPath()+"' was not properly deleted as it should.");
    }

    @Test(groups={"bucket"})
    public void testCopyAndDeleteGCS() throws IOException {
        final String src = publicTestDir + "empty.vcf";
        File dest = createTempFile("copy-empty", ".vcf");
        final String intermediate = BucketUtils.randomRemotePath(getGCPTestStaging(), "test-copy-empty", ".vcf");
        Assert.assertTrue(BucketUtils.isGcsUrl(intermediate), "!BucketUtils.isCloudStorageUrl(intermediate)");
        BucketUtils.copyFile(src, intermediate);
        BucketUtils.copyFile(intermediate, dest.getPath());
        IOUtil.assertFilesEqual(new File(src), dest);
        Assert.assertTrue(BucketUtils.fileExists(intermediate));
        BucketUtils.deleteFile(intermediate);
        Assert.assertFalse(BucketUtils.fileExists(intermediate));
    }

    @Test(groups={"bucket"})
    public void testGetPathOnGcsDirectory() {
        final String dirPath = "gs://bucket/my/dir/";
        final Path pathOnGcs = BucketUtils.getPathOnGcs(dirPath);
        Assert.assertEquals(pathOnGcs.toUri().toString(), dirPath);
    }

    @Test
    public void testCopyAndDeleteHDFS() throws Exception {
        final String src = publicTestDir + "empty.vcf";
        File dest = createTempFile("copy-empty", ".vcf");

        MiniClusterUtils.runOnIsolatedMiniCluster( cluster -> {
            final String intermediate = BucketUtils.randomRemotePath(MiniClusterUtils.getWorkingDir(cluster).toString(), "test-copy-empty", ".vcf");
            Assert.assertTrue(BucketUtils.isHadoopUrl(intermediate), "!BucketUtils.isHadoopUrl(intermediate)");

            BucketUtils.copyFile(src, intermediate);
            BucketUtils.copyFile(intermediate, dest.getPath());
            IOUtil.assertFilesEqual(new File(src), dest);
            Assert.assertTrue(BucketUtils.fileExists(intermediate));
            BucketUtils.deleteFile(intermediate);
            Assert.assertFalse(BucketUtils.fileExists(intermediate));
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

        long fileSize = BucketUtils.fileSize(file1.getAbsolutePath());
        Assert.assertTrue(fileSize > 0);
        long dirSize = BucketUtils.dirSize(new GATKPath(dir.getAbsolutePath()));
        Assert.assertEquals(dirSize, fileSize * 2);
    }

    @Test(groups={"bucket"})
    public void testDirSizeGCS() throws IOException, GeneralSecurityException {
        final String src = publicTestDir + "empty.vcf";
        final String gcsSubDir = BucketUtils.randomRemotePath(getGCPTestStaging(), "dir-", "/");
        final String intermediate = BucketUtils.randomRemotePath(gcsSubDir, "test-copy-empty", ".vcf");
        BucketUtils.copyFile(src, intermediate);
        Assert.assertTrue(BucketUtils.fileExists(intermediate));

        long srcFileSize = BucketUtils.fileSize(src);
        Assert.assertTrue(srcFileSize > 0);
        long intermediateFileSize = BucketUtils.fileSize(intermediate);
        Assert.assertEquals(intermediateFileSize, srcFileSize);
        long intermediateDirSize = BucketUtils.dirSize(new GATKPath(intermediate));
        Assert.assertEquals(intermediateDirSize, srcFileSize);
        long intermediateParentDirSize = BucketUtils.dirSize(new GATKPath(gcsSubDir));
        Assert.assertEquals(intermediateParentDirSize, srcFileSize);

        BucketUtils.deleteFile(intermediate);
        Assert.assertFalse(BucketUtils.fileExists(intermediate));
    }

    @Test
    public void testAddPrefetcher() throws IOException {
        try(final SeekableByteChannel chan = Files.newByteChannel(Paths.get(publicTestDir, "exampleFASTA.fasta"))){
            final SeekableByteChannel prefetchingChannel = BucketUtils.addPrefetcher(1, chan);
            Assert.assertTrue((prefetchingChannel instanceof SeekableByteChannelPrefetcher));
        }
    }

    @DataProvider
    public Object[][] getBufferSizes(){
        return new Object[][]{
                        {0, false},
                        {1, true},
                        {10, true}
                };
    }

    @Test(dataProvider = "getBufferSizes")
    public void testGetWrapper(int bufferSize, boolean prefetchingIsEnabled) throws IOException {
        final Function<SeekableByteChannel, SeekableByteChannel> wrapper = BucketUtils.getPrefetchingWrapper(bufferSize);
        try(final SeekableByteChannel chan = Files.newByteChannel(Paths.get(publicTestDir, "exampleFASTA.fasta"))){
            SeekableByteChannel wrapped = wrapper.apply(chan);
            Assert.assertEquals(wrapped instanceof SeekableByteChannelPrefetcher, prefetchingIsEnabled);
        }
    }

    @Test(groups="cloud")
    public void testCreateSignedUrl() throws IOException {
        final String gcsPathString = getGCPTestInputPath() + "nio/big.txt";
        Assert.assertTrue(Files.exists(IOUtils.getPath(gcsPathString)), "test file is missing, " + gcsPathString);

        final String signed = BucketUtils.createSignedUrlToGcsObject(gcsPathString, 1);
        final Path path = IOUtils.getPath(signed);
        Assert.assertTrue(signed.startsWith("https://"), "path doesn't start with https, " + signed);
        Assert.assertTrue(signed.contains("big.txt"), "path is missing blob name, "+ signed);
        Assert.assertTrue(Files.exists(path), "path doesn't exist: " + signed);
        try(final Stream<String> lines = Files.lines(path))  {
            Assert.assertTrue(lines.anyMatch(line -> line.contains(LINE_IN_REMOTE_TEXT_FILE)), "blob data is incorrect, " + signed);
        }
    }

    @Test(groups="cloud")
    public void testBucketPathToPublicHttpUrl(){
        final String gcsPathString = "gs://hellbender/test/resources/nio/big.txt";
        Assert.assertTrue(Files.exists(IOUtils.getPath(gcsPathString)), "Test file is missing: "+ gcsPathString);
        final String publicHttpLink = BucketUtils.bucketPathToPublicHttpUrl(gcsPathString);
        Assert.assertEquals(publicHttpLink, "https://storage.googleapis.com/hellbender/test/resources/nio/big.txt");
        Assert.assertTrue(Files.exists(IOUtils.getPath(publicHttpLink)));
    }

    @Test(groups="cloud", singleThreaded=true)
    public void testRequesterPays() throws IOException {

        // IOUtils.getPath() triggers the requester pays check, it's not lazily delayed until trying to read
        // the file.
        try {
            //Assert that this fails when no project is provided
            setNoProjectForRequesterPays();
            Assert.assertThrows(StorageException.class, () -> {
                try( final Stream<String> lines = Files.lines(IOUtils.getPath(FILE_IN_REQUESTER_PAYS_BUCKET))) {
                   lines.anyMatch(line -> line.contains(LINE_IN_REMOTE_TEXT_FILE));
                }
            });

            setRequesterPays();
            try( final Stream<String> lines = Files.lines(IOUtils.getPath(FILE_IN_REQUESTER_PAYS_BUCKET))) {
                Assert.assertTrue(lines.anyMatch(line -> line.contains(LINE_IN_REMOTE_TEXT_FILE)),
                        "Failed to read from file.");
            }
        } finally {
            setDefaultNioOptions();
        }
    }

}
