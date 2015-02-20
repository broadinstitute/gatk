/*
* Copyright (c) 2012 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.io;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.Random;

public class IOUtilsUnitTest extends BaseTest {


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
        File nonExistentFile = new File("djfhsdkjghdfk");
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
}
