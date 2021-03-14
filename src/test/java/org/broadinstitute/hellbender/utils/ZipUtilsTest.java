package org.broadinstitute.hellbender.utils;

import com.google.common.io.Files;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Tests cases for {@link ZipUtils}.
 */
public final class ZipUtilsTest {

    private static int AVERAGE_FILE_SIZE_IN_BYTES = 1024;
    private static double SD_FILE_SIZE = 0.5;

    @Test()
    public void testZipAndUnzipSingleFile() throws IOException {
        final File root = createCase(new String[] {"a.bin"}, new Random(13));
        final File aBin = new File(root, "a.bin");
        Assert.assertTrue(aBin.isFile());
        final File zipFile = File.createTempFile("test", ".zip");
        zipFile.delete();
        Assert.assertFalse(zipFile.exists());
        ZipUtils.zip(aBin, new GATKPath(zipFile.toString()));
        Assert.assertTrue(zipFile.isFile());
        final File unzipRoot = Files.createTempDir();
        unzipRoot.delete();
        Assert.assertFalse(unzipRoot.exists());
        ZipUtils.unzip(new GATKPath(zipFile.toString()), unzipRoot);
        Assert.assertTrue(unzipRoot.isDirectory());
        final File[] subFiles = unzipRoot.listFiles();
        Assert.assertEquals(subFiles.length, 1);
        Assert.assertEquals(subFiles[0].toString().replaceAll(unzipRoot.toString() + "/?", ""), "a.bin");
        assertEqualContents(aBin, subFiles[0]);
        FileUtils.deleteDirectory(root);
    }

    @Test(dataProvider="zipAndUnzipFullFileData")
    public void testZipAndUnzipFullFile(final File root, final String[] files) throws IOException {
        final File zipFile = File.createTempFile("test", ".zip");
        Assert.assertTrue(zipFile.delete());
        Assert.assertFalse(zipFile.exists());
        ZipUtils.zip(root, new GATKPath(zipFile.toString()));
        zipFile.deleteOnExit();
        Assert.assertTrue(zipFile.exists());
        Assert.assertTrue(zipFile.isFile());
        // now unzip.
        final File unzipRoot = Files.createTempDir();
        FileUtils.deleteDirectory(unzipRoot);
        Assert.assertFalse(unzipRoot.exists());
        ZipUtils.unzip(new GATKPath(zipFile.toString()), unzipRoot);
        Assert.assertTrue(unzipRoot.exists());
        // now we compare the unzipped content with the original input content.
        for (final String inFileName : files) {
            final File inFile = new File(root, inFileName);
            final File outFile = new File(unzipRoot, inFileName);
            assertEqualContents(inFile, outFile);
        }
    }

    @Test(dataProvider="zipAndUnzipFullFileData")
    public void testZipFullAndUnzipSubset(final File root, final String[] files) throws IOException {
        final File zipFile = File.createTempFile("test", ".zip");
        zipFile.delete();
        zipFile.deleteOnExit();
        Assert.assertFalse(zipFile.exists());
        ZipUtils.zip(root, new GATKPath(zipFile.toString()));
        Assert.assertTrue(zipFile.exists());
        Assert.assertTrue(zipFile.isFile());
        // now unzip.
        final File unzipRoot = Files.createTempDir();
        FileUtils.deleteDirectory(unzipRoot);
        Assert.assertFalse(unzipRoot.exists());
        final Random rdn = new Random(Arrays.hashCode(files));
        final String[] only = IntStream.range(0, files.length)
           .filter(i -> rdn.nextBoolean())
           .mapToObj(i -> files[i])
           .toArray(String[]::new);

        ZipUtils.unzip(new GATKPath(zipFile.toString()), unzipRoot, only);
        Assert.assertTrue(unzipRoot.exists());
        // now we compare the unzipped content with the original input content.
        for (final String inFileName : files) {
            final File inFile = new File(root, inFileName);
            final File outFile = new File(unzipRoot, inFileName);
            if (only.length == 0 || Arrays.stream(only).anyMatch(s -> s.equals(inFileName))) {
                assertEqualContents(inFile, outFile);
            } else {
                Assert.assertFalse(outFile.exists(), outFile.toString());
            }
        }
    }

    @Test(dataProvider="zipAndUnzipFullFileData")
    public void testZipSubsetAndUnzipFull(final File root, final String[] files) throws IOException {
        final File zipFile = File.createTempFile("test", ".zip");
        zipFile.delete();
        zipFile.deleteOnExit();
        Assert.assertFalse(zipFile.exists());
        final Random rdn = new Random(Arrays.hashCode(files) * 31);
        final String[] only = IntStream.range(0, files.length)
                .filter(i -> rdn.nextBoolean())
                .mapToObj(i -> files[i])
                .toArray(String[]::new);

        ZipUtils.zip(root, new GATKPath(zipFile.toString()), only);
        Assert.assertTrue(zipFile.exists());
        Assert.assertTrue(zipFile.isFile());
        // now unzip.
        final File unzipRoot = Files.createTempDir();
        FileUtils.deleteDirectory(unzipRoot);
        Assert.assertFalse(unzipRoot.exists());

        ZipUtils.unzip(new GATKPath(zipFile.toString()), unzipRoot);
        Assert.assertTrue(unzipRoot.exists());
        // now we compare the unzipped content with the original input content.
        for (final String inFileName : files) {
            final File inFile = new File(root, inFileName);
            final File outFile = new File(unzipRoot, inFileName);
            if (only.length == 0 || Arrays.stream(only).anyMatch(s -> s.equals(inFile.toString()))) {
                assertEqualContents(inFile, outFile);
            } else {
                Assert.assertFalse(outFile.exists(), outFile.toString());
            }
        }
    }

    @DataProvider
    public Object[][] zipAndUnzipFullFileData() {
        final Random rdn = new Random(13);
        final String[] cases = {
                "a/b/c.txt,a/b/d.txt,a/b.d,a/b.e/,b,c,d,e,e/a.txt",
                "a,b,c,d/,e/a,e/b",
                "a.txt",
                "a/,a/a.txt/",
        };
        return Arrays.stream(cases)
                .map(s -> s.split(","))
                .map(ss -> new Object[] { createCase(ss, rdn), ss })
                .toArray(Object[][]::new);
    }


    private void assertEqualContents(final File inFile, final File outFile) throws IOException {
        Assert.assertTrue(inFile.exists()); // paranoia check, it should be true.
        Assert.assertTrue(outFile.exists(), outFile.toString());
        Assert.assertEquals(inFile.isFile(), outFile.isFile());
        Assert.assertEquals(inFile.isDirectory(), outFile.isDirectory());
        Assert.assertEquals(inFile.length(), outFile.length());
        if (inFile.isDirectory()) {
            final File[] inSubFiles = Arrays.asList(inFile.listFiles()).stream()
                    .map(f -> new File(f.toString().replaceAll(inFile.toString() + "/?", "")))
                    .toArray(File[]::new);
            final File[] outSubFiles = Arrays.asList(outFile.listFiles()).stream()
                    .map(f -> new File(f.toString().replaceAll(outFile.toString() + "/?", "")))
                    .toArray(File[]::new);
            Arrays.sort(inSubFiles);
            Arrays.sort(outSubFiles);
            Assert.assertEquals(inSubFiles, outSubFiles);
        } else {
            final byte[] inBuffer = fileIntoBytes(inFile);
            final byte[] outBuffer = fileIntoBytes(outFile);
            Assert.assertEquals(inBuffer, outBuffer);
        }

    }

    private byte[] fileIntoBytes(final File file) throws IOException {
        final int length = (int) file.length();
        final byte[] result = new byte[length];
        try (final InputStream is = new FileInputStream(file)) {
           is.read(result);
        }
        return result;
    }

    private File createCase(final String[] files, final Random rdn) throws UncheckedIOException {
        try {
            final File root = File.createTempFile("dir", ".dir");
            root.delete();
            root.mkdir();
            final List<String> fileNames = new ArrayList<>(Arrays.asList(files));
            final Comparator<String> comparator = Comparator.reverseOrder();
            fileNames.sort(comparator);
            for (final String fileName : fileNames) {
                final File file = new File(root, fileName);
                if (fileName.endsWith(File.pathSeparator)) {
                    file.mkdirs();
                } else if (!file.exists()) {
                    file.getParentFile().mkdirs();
                    file.createNewFile();
                    fillFileWithGarbage(file, rdn);
                }
            }
            return root;
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }

    private void fillFileWithGarbage(final File dest, final Random rdn) throws IOException {
        final int fileSize = Math.max(100, (int) Math.round(AVERAGE_FILE_SIZE_IN_BYTES + rdn.nextGaussian() * AVERAGE_FILE_SIZE_IN_BYTES * SD_FILE_SIZE));
        try (final OutputStream os = new FileOutputStream(dest)) {
            if (fileSize < 0) {
                throw new RuntimeException("blah");
            }
            final byte[] content = new byte[fileSize];
            rdn.nextBytes(content);
            os.write(content);
        }
    }
}
