package org.broadinstitute.hellbender.utils.io;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import org.apache.logging.log4j.core.util.FileUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.*;
import java.util.Arrays;
import java.util.Random;
import java.util.UUID;

public final class IOUtilsUnitTest extends GATKBaseTest {

    private FileSystem jimfs;

    @BeforeClass
    public void setUp() {
        jimfs = Jimfs.newFileSystem(Configuration.unix());
    }
    @AfterClass
    public void tearDown() {
        try {
            jimfs.close();
        } catch (final IOException e) {
            log("Unable to close JimFS");
        }
    }

    @Test
    public void testTempDir() {
        File tempDir = IOUtils.createTempDir("Q-Unit-Test");
        Assert.assertTrue(tempDir.exists());
        Assert.assertFalse(tempDir.isFile());
        Assert.assertTrue(tempDir.isDirectory());
        boolean deleted = IOUtils.tryDelete(tempDir);
        Assert.assertTrue(deleted);
        Assert.assertFalse(tempDir.exists());
    }

    @Test
    public void testIsSpecialFile() {
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/null")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/full")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/stdout")));
        Assert.assertTrue(IOUtils.isSpecialFile(new File("/dev/stderr")));
        Assert.assertFalse(IOUtils.isSpecialFile(null));
        Assert.assertFalse(IOUtils.isSpecialFile(new File("/home/user/my.file")));
        Assert.assertFalse(IOUtils.isSpecialFile(new File("/devfake/null")));
    }

    @DataProvider( name = "ByteArrayIOTestData")
    public Object[][] byteArrayIOTestDataProvider() {
        return new Object[][] {
                // file size, read buffer size
                { 0,     4096 },
                { 1,     4096 },
                { 2000,  4096 },
                { 4095,  4096 },
                { 4096,  4096 },
                { 4097,  4096 },
                { 6000,  4096 },
                { 8191,  4096 },
                { 8192,  4096 },
                { 8193,  4096 },
                { 10000, 4096 }
        };
    }

    @Test( dataProvider = "ByteArrayIOTestData" )
    public void testWriteThenReadFileIntoByteArray ( int fileSize, int readBufferSize ) throws Exception {
        File tempFile = createTempFile(String.format("testWriteThenReadFileIntoByteArray_%d_%d", fileSize, readBufferSize), "tmp");

        byte[] dataWritten = getDeterministicRandomData(fileSize);
        IOUtils.writeByteArrayToFile(dataWritten, tempFile);
        byte[] dataRead = IOUtils.readFileIntoByteArray(tempFile, readBufferSize);

        Assert.assertEquals(dataRead.length, dataWritten.length);
        Assert.assertTrue(Arrays.equals(dataRead, dataWritten));
    }

    @Test( dataProvider = "ByteArrayIOTestData" )
    public void testWriteThenReadStreamIntoByteArray ( int fileSize, int readBufferSize ) throws Exception {
        File tempFile = createTempFile(String.format("testWriteThenReadStreamIntoByteArray_%d_%d", fileSize, readBufferSize), "tmp");

        byte[] dataWritten = getDeterministicRandomData(fileSize);
        IOUtils.writeByteArrayToStream(dataWritten, new FileOutputStream(tempFile));
        byte[] dataRead = IOUtils.readStreamIntoByteArray(new FileInputStream(tempFile), readBufferSize);

        Assert.assertEquals(dataRead.length, dataWritten.length);
        Assert.assertTrue(Arrays.equals(dataRead, dataWritten));
    }

    @Test( expectedExceptions = UserException.CouldNotReadInputFile.class )
    public void testReadNonExistentFileIntoByteArray() {
        File nonExistentFile = GATKBaseTest.getSafeNonExistentFile("djfhsdkjghdfk");
        Assert.assertFalse(nonExistentFile.exists());

        IOUtils.readFileIntoByteArray(nonExistentFile);
    }

    @Test( expectedExceptions = IllegalArgumentException.class )
    public void testReadStreamIntoByteArrayInvalidBufferSize() throws Exception {
        IOUtils.readStreamIntoByteArray(new FileInputStream(createTempFile("testReadStreamIntoByteArrayInvalidBufferSize", "tmp")),
                -1);
    }

    private byte[] getDeterministicRandomData ( int size ) {
        Utils.resetRandomGenerator();
        Random rand = Utils.getRandomGenerator();

        byte[] randomData = new byte[size];
        rand.nextBytes(randomData);

        return randomData;
    }

    @Test
    public void testDeleteDirOnExit() throws IOException {
        //This just tests that the code runs without crashing.
        //It runs at jvm shutdown so there isn't a good way to test it properly.
        //If you see a directory in the hellbender main folder called

        final File dir = new File(GATKBaseTest.publicTestDir + "I_SHOULD_HAVE_BEEN_DELETED");
        IOUtils.deleteRecursivelyOnExit(dir);

        FileUtils.mkdir(dir, true);
        final File subdir = new File(dir, "subdir");
        FileUtils.mkdir(subdir, true);
        File someFile = new File(dir, "someFile");
        someFile.createNewFile();
        File anotherFile = new File(subdir, "anotherFile");
        anotherFile.createNewFile();
    }

    @DataProvider(name = "extensionsToReplace")
    public Object[][] getExtensionsToReplace(){
        return new Object[][] {
                {"file.old", "file.new"},
                {"file.something.old", "file.something.new"},
                {"src/test.something/file", "src/test.something/file.new"},
                {"/.src.folder/some/thing/.secret/file.old", "/.src.folder/some/thing/.secret/file.new" }
        };
    }

    @Test(dataProvider = "extensionsToReplace")
    public void testReplaceExtension(String input, String expected){
        Assert.assertEquals(IOUtils.replaceExtension(input, "new"), expected);
        Assert.assertEquals(IOUtils.replaceExtension(input, "new"), IOUtils.replaceExtension(input,"..new"));
        Assert.assertEquals(IOUtils.replaceExtension(new File(input), "new"), new File(expected));
    }

    @Test(groups={"bucket"})
    public void testGetPath() throws IOException {
        innerTestGetPath(getGCPTestInputPath() + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");
        innerTestGetPath("file://" + NA12878_20_21_WGS_bam);
        innerTestGetPath(NA12878_20_21_WGS_bam);
    }

    @Test
    public void testGetPathHandlesIntervals() throws IOException {
        // Make sure we don't crash if passing intervals to getPath.
        // Also, it shouldn't crash when we check whether the file exists.
        Assert.assertFalse(Files.exists(IOUtils.getPath("chr1:10-11")));
        Assert.assertFalse(Files.exists(IOUtils.getPath("chr1:10")));
        Assert.assertFalse(Files.exists(IOUtils.getPath("1:10-11")));
        Assert.assertFalse(Files.exists(IOUtils.getPath("1:10")));
    }

    private void innerTestGetPath(String s) throws IOException {
        Path p = IOUtils.getPath(s);
        long size = Files.size(p);
        Assert.assertTrue(size>0);
    }

    @DataProvider
    public Object[][] absoluteNames() {
        return new Object[][] {
                {"relative/example.txt", new File("relative/example.txt").getAbsolutePath()},
                {"/local/example.txt", "/local/example.txt"},
                // note that path normalization removes the extra / in file://
                {"/local/file://example.txt", "/local/file:/example.txt"},
                {"file:///local/example.txt", "/local/example.txt"},
                {"gs://dir/example.txt", "gs://dir/example.txt"}
        };
    }

    @Test(dataProvider = "absoluteNames")
    public void testGetAbsolutePathWithoutFileProtocol(final String uriString, final String expected) {
        Assert.assertEquals(IOUtils.getAbsolutePathWithoutFileProtocol(IOUtils.getPath(uriString)), expected);
    }

    @Test
    public void testAppendPathToDir() throws Exception {
        Assert.assertEquals(IOUtils.appendPathToDir("dir", "file"), "dir/file");
        Assert.assertEquals(IOUtils.appendPathToDir("dir/", "file"), "dir/file");
        Assert.assertEquals(IOUtils.appendPathToDir("dir", "/file"), "/file");
        Assert.assertEquals(IOUtils.appendPathToDir("dir/", "/file"), "/file");
        Assert.assertEquals(IOUtils.appendPathToDir("/path/to/dir", "anotherdir/file"), "/path/to/dir/anotherdir/file");

        // hdfs: URI
        Path tempPath = IOUtils.getPath(MiniClusterUtils.getWorkingDir(MiniClusterUtils.getMiniCluster()).toUri().toString());
        Assert.assertEquals(IOUtils.appendPathToDir(tempPath.toString(), "temp"), tempPath.toString()+"/temp");

        // gs: URI
        Assert.assertEquals(IOUtils.appendPathToDir("gs://abucket/dir", "file"), "gs://abucket/dir/file");

        // file: URI
        Assert.assertEquals(IOUtils.appendPathToDir("file:///dir", "file"), "file:///dir/file");
    }

    @Test
    public void testSuccessfulCanReadFileCheck() {
        final File expectedFile = createTempFile("Utils-can-read-test",".txt");
        IOUtils.canReadFile(expectedFile);
    }

    @Test
    public void testSuccessfulCanReadFilesCheck() {
        final File file1 = createTempFile("Utils-can-read-test1",".txt");
        final File file2 = createTempFile("Utils-can-read-test2",".txt");
        IOUtils.canReadFile(file1, file2);
    }

    @DataProvider(name = "unsuccessfulCanReadFileCheckData")
    public Object[][] unsuccessfulCanReadFileCheckData() {
        final File directory = createTempDir("Utils-can-read-file-Dir");
        final File nonExistingFile = createTempFile("Utils-cant-read-NoFile", ".file");
        nonExistingFile.delete();
        final File nonReadable = createTempFile("Utils-cant-read-NotReadable", ".file");

        //Note: if this test suite is run as root (eg in a Docker image), setting readable to false may fail.
        // So we check it here and skip this test if we can't make a file non readable.
        final boolean successSetNotReadable = nonReadable.setReadable(false) && !nonReadable.canRead();
        if (successSetNotReadable) {
            return new Object[][]{{directory}, {nonExistingFile}, {nonReadable}};
        } else {
            //pass in a special value to be able to throw a SkipException and make the same number of tests
            logger.debug("cannot make a file unreadable (maybe you're running as root)");
            return new Object[][]{{directory}, {nonExistingFile}, {null}};
        }
    }

    @Test(dataProvider = "unsuccessfulCanReadFileCheckData",
            expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testUnsuccessfulCanReadFileCheck(final File file) {
        if (file == null){
            throw new SkipException("cannot make a file unreadable (maybe you're running as root)");
        }
        IOUtils.canReadFile(file);
    }

    @Test
    public void testCreateDirectory() throws IOException {

        // hdfs
        Path tempPath = IOUtils.getPath(MiniClusterUtils.getWorkingDir(MiniClusterUtils.getMiniCluster()).toUri().toString())
                .resolve("temp");
        IOUtils.createDirectory(tempPath.toUri().toString());
        Assert.assertTrue(java.nio.file.Files.exists(tempPath));

        // local
        tempPath = IOUtils.getPath(BaseTest.createTempDir("xxx").getAbsolutePath().concat("/sub"));
        IOUtils.createDirectory(tempPath.toUri().toString());
        Assert.assertTrue(java.nio.file.Files.exists(tempPath));
    }

    @DataProvider(name="urlEncodeDecode")
    public Object[][] urlEncodeDecode() {
        return new Object[][]{
                // string, url encoding
                { "string", "string"},
                { "string.", "string."},
                { "string1", "string1"},
                { "string with space", "string+with+space"},
                { "string://", "string%3A%2F%2F"},
        };
    }

    @Test(dataProvider = "urlEncodeDecode")
    public void testUrlEncodeDecode(final String string, final String encoded) {
        Assert.assertEquals(IOUtils.urlEncode(string), encoded);
        Assert.assertEquals(string, IOUtils.urlDecode(encoded));
    }

    @DataProvider(name="resourcePaths")
    public Object[][] getResourcePaths() {
        final String testResourceContents = "this is a test resource file";
        return new Object[][] {
                { Resource.LARGE_RUNTIME_RESOURCES_PATH + "/testResourceFile.txt", null , testResourceContents},
                { "testResourceFile.txt", IOUtils.class , testResourceContents},
        };
    }

    @Test(dataProvider = "resourcePaths")
    public void testWriteTempResource(
            final String resourcePath, final Class<?> relativeClass, final String expectedFirstLine) throws IOException
    {
        final Resource largeResource = new Resource(resourcePath, relativeClass);
        final File resourceFile = IOUtils.writeTempResource(largeResource);
        final String resourceContentsFirstLine = getFirstLineAndDeleteTempFile(resourceFile);
        Assert.assertEquals(resourceContentsFirstLine, expectedFirstLine);
    }

    @Test(dataProvider = "resourcePaths")
    public void testWriteTempResourceFromPath(
            final String resourcePath, final Class<?> relativeClass, final String expectedFirstLine) throws IOException
    {
        final File resourceFile = IOUtils.writeTempResourceFromPath(resourcePath, relativeClass);
        final String resourceContentsFirstLine = getFirstLineAndDeleteTempFile(resourceFile);
        Assert.assertEquals(resourceContentsFirstLine, expectedFirstLine);
    }

    @Test
    public void testCreateTempFileInDirectory() {
        final File tempDir = createTempDir("testCreateTempFileInDirectory");
        final File tempFile = IOUtils.createTempFileInDirectory("testCreateTempFileInDirectory", ".txt", tempDir);
        Assert.assertTrue(tempFile.exists(), "file was not written to temp file: " + tempFile);
        Assert.assertEquals(tempFile.getParentFile().getAbsolutePath(), (tempDir.getAbsolutePath()),
                "file was not written to temp file: " + tempFile + " in dir: " + tempDir);
    }

    @DataProvider
    public Object[][] tmpPathDirs() throws Exception {
        return new Object[][] {
                {createTempDir("local").toPath()},
                {Files.createDirectory(jimfs.getPath("tmp"))}
        };
    }

    @Test(dataProvider = "tmpPathDirs", singleThreaded = true)
    public void testTempPath(final Path tempDir) throws Exception {
        // store the previous tmp.dir to check that it is working
        final String previousTmpDir = System.getProperty("java.io.tmpdir");
        try {
            System.setProperty("java.io.tmpdir", IOUtils.getAbsolutePathWithoutFileProtocol(tempDir));
            final Path tempFile = IOUtils.createTempPath(UUID.randomUUID().toString(), ".txt");
            Assert.assertTrue(Files.exists(tempFile),
                    "file was not written to temp file: " + tempFile);
            Assert.assertEquals(tempFile.getParent().toUri().toString(), tempDir.toUri().toString(),
                    "file was not written to temp file: " + tempFile + " in dir: " + tempDir);
        } finally {
            System.setProperty("java.io.tmpdir", previousTmpDir);
        }
    }

    private String getFirstLineAndDeleteTempFile(final File tempResourceFile) throws IOException {
        String resourceContentsFirstLine = "";
        try (final FileReader fr = new FileReader(tempResourceFile);
             final BufferedReader br = new BufferedReader(fr)) {
            resourceContentsFirstLine = br.readLine();
        }
        finally {
            tempResourceFile.delete();
        }
        return resourceContentsFirstLine;
    }

    @DataProvider
    private Object[][] getIsHDF5TestFiles() {
        return new Object[][] {
                { getToolTestDataDir() + "/isValidHDF5.hdf5", true },
                { getToolTestDataDir() + "isValidHDF5.ext", true },
                { getToolTestDataDir() + "/isTSV.tsv", false },
                { getToolTestDataDir() + "/isTSV.hdf", false },
                { getToolTestDataDir() + "/isTSV.hdf5", false },
                { getToolTestDataDir() + "/isEmpty.txt", false },
        };
    }

    @Test(dataProvider = "getIsHDF5TestFiles")
    public void testIsHDF5File(final String filePath, final boolean expected) {
        final Path testPath = Paths.get(filePath);
        Assert.assertEquals(IOUtils.isHDF5File(testPath), expected);
    }

}
