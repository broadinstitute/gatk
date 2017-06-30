package org.broadinstitute.hellbender.utils.io;

import org.apache.logging.log4j.core.util.FileUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Random;

public final class IOUtilsUnitTest extends BaseTest {

    @Test
    public void testTempDir() {
        File tempDir = IOUtils.tempDir("Q-Unit-Test", "", new File("queueTempDirToDelete"));
        Assert.assertTrue(tempDir.exists());
        Assert.assertFalse(tempDir.isFile());
        Assert.assertTrue(tempDir.isDirectory());
        boolean deleted = IOUtils.tryDelete(tempDir);
        Assert.assertTrue(deleted);
        Assert.assertFalse(tempDir.exists());
    }

    @Test
    public void testAbsolute() {
        File dir = IOUtils.absolute(new File("/path/./to/./directory/."));
        Assert.assertEquals(dir, new File("/path/to/directory"));

        dir = IOUtils.absolute(new File("/"));
        Assert.assertEquals(dir, new File("/"));

        dir = IOUtils.absolute(new File("/."));
        Assert.assertEquals(dir, new File("/"));

        dir = IOUtils.absolute(new File("/././."));
        Assert.assertEquals(dir, new File("/"));

        dir = IOUtils.absolute(new File("/./directory/."));
        Assert.assertEquals(dir, new File("/directory"));

        dir = IOUtils.absolute(new File("/./directory/./"));
        Assert.assertEquals(dir, new File("/directory"));

        dir = IOUtils.absolute(new File("/./directory./"));
        Assert.assertEquals(dir, new File("/directory."));

        dir = IOUtils.absolute(new File("/./.directory/"));
        Assert.assertEquals(dir, new File("/.directory"));
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
        File nonExistentFile = BaseTest.getSafeNonExistentFile("djfhsdkjghdfk");
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

        final File dir = new File(TestResources.publicTestDir + "I_SHOULD_HAVE_BEEN_DELETED");
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
        innerTestGetPath("file://" + TestResources.NA12878_20_21_WGS_bam);
        innerTestGetPath(TestResources.NA12878_20_21_WGS_bam);
    }

    private void innerTestGetPath(String s) throws IOException {
        Path p = IOUtils.getPath(s);
        long size = Files.size(p);
        Assert.assertTrue(size>0);
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
}
