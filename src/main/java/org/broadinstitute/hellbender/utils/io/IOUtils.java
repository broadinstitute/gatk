package org.broadinstitute.hellbender.utils.io;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.util.TabixUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.net.URI;
import java.nio.file.*;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

public final class IOUtils {
    private static final Logger logger = LogManager.getLogger(IOUtils.class);
    private static final File DEV_DIR = new File("/dev");

    /**
     * Returns true if the file's extension is CRAM.
     */
    public static boolean isCramFile(final File inputFile) {
        return isCramFileName(inputFile.getName());
    }

    /**
     * Returns true if the file's extension is CRAM.
     */
    public static boolean isCramFileName(final String inputFileName) {
        return CramIO.CRAM_FILE_EXTENSION.equalsIgnoreCase("." + FilenameUtils.getExtension(inputFileName));
    }

    /**
     * Returns true if the file's extension is BAM.
     */
    public static boolean isBamFileName(final String inputFileName) {
        return BamFileIoUtils.BAM_FILE_EXTENSION.equalsIgnoreCase("." + FilenameUtils.getExtension(inputFileName));
    }

    /**
     * Creates a temp directory with the prefix and optional suffix.
     *
     * @param prefix       Prefix for the directory name.
     * @param suffix       Optional suffix for the directory name.
     * @return The created temporary directory.
     */
    public static File tempDir(String prefix, String suffix) {
        return tempDir(prefix, suffix, null);
    }

    /**
     * Creates a temp directory with the prefix and optional suffix.
     *
     * @param prefix        Prefix for the directory name.
     * @param suffix        Optional suffix for the directory name.
     * @param tempDirParent Parent directory for the temp directory.
     * @return The created temporary directory.
     */
    public static File tempDir(String prefix, String suffix, File tempDirParent) {
        try {
            if (tempDirParent == null)
                tempDirParent = FileUtils.getTempDirectory();
            if (!tempDirParent.exists() && !tempDirParent.mkdirs())
                throw new UserException.BadTmpDir("Could not create temp directory: " + tempDirParent);
            File temp = File.createTempFile(prefix, suffix, tempDirParent);
            if (!temp.delete())
                throw new UserException.BadTmpDir("Could not delete sub file: " + temp.getAbsolutePath());
            if (!temp.mkdir())
                throw new UserException.BadTmpDir("Could not create sub directory: " + temp.getAbsolutePath());
            return absolute(temp);
        } catch (IOException e) {
            throw new UserException.BadTmpDir(e.getMessage());
        }
    }

    /**
     * Writes content to a temp file and returns the path to the temporary file.
     *
     * @param content   to write.
     * @param prefix    Prefix for the temp file.
     * @param suffix    Suffix for the temp file.
     * @return the path to the temp file.
     */
    public static File writeTempFile(String content, String prefix, String suffix) {
        return writeTempFile(content, prefix, suffix, null);
    }

    /**
     * Writes content to a temp file and returns the path to the temporary file.
     *
     * @param content   to write.
     * @param prefix    Prefix for the temp file.
     * @param suffix    Suffix for the temp file.
     * @param directory Directory for the temp file.
     * @return the path to the temp file.
     */
    public static File writeTempFile(String content, String prefix, String suffix, File directory) {
        try {
            File tempFile = absolute(File.createTempFile(prefix, suffix, directory));
            FileUtils.writeStringToFile(tempFile, content);
            return tempFile;
        } catch (IOException e) {
            throw new UserException.BadTmpDir(e.getMessage());
        }
    }

    /**
     * Returns true if the file is a special file.
     * @param file File path to check.
     * @return true if the file is a special file.
     */
    public static boolean isSpecialFile(File file) {
        return file != null && (file.getAbsolutePath().startsWith("/dev/") || file.equals(DEV_DIR));
    }

    /**
     * Tries to delete a file. Emits a warning if the file
     * is not a special file and was unable to be deleted.
     *
     * @param file File to delete.
     * @return true if the file was deleted.
     */
    public static boolean tryDelete(File file) {
        if (isSpecialFile(file)) {
            logger.debug("Not trying to delete " + file);
            return false;
        }
        boolean deleted = FileUtils.deleteQuietly(file);
        if (deleted)
            logger.debug("Deleted " + file);
        else if (file.exists())
            logger.warn("Unable to delete " + file);
        return deleted;
    }

    /**
     * A mix of getCanonicalFile and getAbsoluteFile that returns the
     * absolute path to the file without deferencing symbolic links.
     *
     * @param file the file.
     * @return the absolute path to the file.
     */
    public static File absolute(File file) {
        return replacePath(file, absolutePath(file));
    }

    private static String absolutePath(File file) {
        File fileAbs = file.getAbsoluteFile();
        LinkedList<String> names = new LinkedList<>();
        while (fileAbs != null) {
            String name = fileAbs.getName();
            fileAbs = fileAbs.getParentFile();

            if (".".equals(name)) {
                /* skip */

                /* TODO: What do we do for ".."?
              } else if (name == "..") {

                CentOS tcsh says use getCanonicalFile:
                ~ $ mkdir -p test1/test2
                ~ $ ln -s test1/test2 test3
                ~ $ cd test3/..
                ~/test1 $

                Mac bash says keep going with getAbsoluteFile:
                ~ $ mkdir -p test1/test2
                ~ $ ln -s test1/test2 test3
                ~ $ cd test3/..
                ~ $

                For now, leave it and let the shell figure it out.
                */
            } else {
                names.add(0, name);
            }
        }

        return ("/" + StringUtils.join(names, "/"));
    }

    private static File replacePath(File file, String path) {
        if (file instanceof FileExtension)
            return ((FileExtension)file).withPath(path);
        if (!File.class.equals(file.getClass()))
            throw new GATKException("Sub classes of java.io.File must also implement FileExtension");
        return new File(path);
    }

    /**
     * Writes the an embedded resource to a temp file.
     * File is not scheduled for deletion and must be cleaned up by the caller.
     * @param resource Embedded resource.
     * @return Path to the temp file with the contents of the resource.
     */
    public static File writeTempResource(Resource resource) {
        File temp;
        try {
            temp = File.createTempFile(FilenameUtils.getBaseName(resource.getPath()) + ".", "." + FilenameUtils.getExtension(resource.getPath()));
        } catch (IOException e) {
            throw new UserException.BadTmpDir(e.getMessage());
        }
        writeResource(resource, temp);
        return temp;
    }

    /**
     * Writes the an embedded resource to a file.
     * File is not scheduled for deletion and must be cleaned up by the caller.
     * @param resource Embedded resource.
     * @param file File path to write.
     */
    public static void writeResource(Resource resource, File file) {
        String path = resource.getPath();
        InputStream inputStream = resource.getResourceContentsAsStream();
        OutputStream outputStream = null;
        try {
            outputStream = FileUtils.openOutputStream(file);
            org.apache.commons.io.IOUtils.copy(inputStream, outputStream);
        } catch (IOException e) {
            throw new GATKException(String.format("Unable to copy resource '%s' to '%s'", path, file), e);
        } finally {
            org.apache.commons.io.IOUtils.closeQuietly(inputStream);
            org.apache.commons.io.IOUtils.closeQuietly(outputStream);
        }
    }

    /**
     * Reads the entirety of the given file into a byte array. Uses a read buffer size of 4096 bytes.
     *
     * @param source File to read
     * @return The contents of the file as a byte array
     */
    public static byte[] readFileIntoByteArray ( File source ) {
        return readFileIntoByteArray(source, 4096);
    }

    /**
     * Reads the entirety of the given file into a byte array using the requested read buffer size.
     *
     * @param source File to read
     * @param readBufferSize Number of bytes to read in at one time
     * @return The contents of the file as a byte array
     */
    public static byte[] readFileIntoByteArray ( File source, int readBufferSize ) {
        if ( source == null ) {
            throw new GATKException("Source file was null");
        }

        byte[] fileContents;

        try {
            fileContents = readStreamIntoByteArray(new FileInputStream(source), readBufferSize);
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(source, e);
        }

        if ( fileContents.length != source.length() ) {
            throw new UserException.CouldNotReadInputFile(String.format("Unable to completely read file %s: read only %d/%d bytes",
                    source.getAbsolutePath(), fileContents.length, source.length()));
        }

        return fileContents;
    }

    /**
     * Reads all data from the given stream into a byte array using the requested read buffer size.
     *
     * @param in Stream to read data from
     * @param readBufferSize Number of bytes to read in at one time
     * @return The contents of the stream as a byte array
     */
    public static byte[] readStreamIntoByteArray ( InputStream in, int readBufferSize ) {
        if ( in == null ) {
            throw new IllegalArgumentException("Input stream was null");
        }
        else if ( readBufferSize <= 0 ) {
            throw new IllegalArgumentException("Read buffer size must be > 0");
        }

        // Use a fixed-size buffer for each read, but a dynamically-growing buffer
        // to hold the accumulated contents of the file/stream:
        byte[] readBuffer = new byte[readBufferSize];
        ByteArrayOutputStream fileBuffer = new ByteArrayOutputStream(readBufferSize * 4);

        try {
            try {
                int currentBytesRead;

                while ( (currentBytesRead = in.read(readBuffer, 0, readBuffer.length)) >= 0 ) {
                    fileBuffer.write(readBuffer, 0, currentBytesRead);
                }
            }
            finally {
                in.close();
            }
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error reading from input stream", e);
        }

        return fileBuffer.toByteArray();
    }

    /**
     * Writes the given array of bytes to a file
     *
     * @param bytes Data to write
     * @param destination File to write the data to
     */
    public static void writeByteArrayToFile ( byte[] bytes, File destination ) {
        if ( destination == null ) {
            throw new GATKException("Destination file was null");
        }

        try {
            writeByteArrayToStream(bytes, new FileOutputStream(destination));
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(destination, e);
        }
    }

    /**
     * Writes the given array of bytes to a stream
     *
     * @param bytes Data to write
     * @param out Stream to write the data to
     */
    public static void writeByteArrayToStream ( byte[] bytes, OutputStream out ) {
        if ( bytes == null || out == null ) {
            throw new GATKException("Data to write or output stream was null");
        }

        try {
            try {
                out.write(bytes);
            }
            finally {
                out.close();
            }
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile("I/O error writing to output stream", e);
        }
    }

    /**
     * Un-gzips the input file to the output file.
     */
    public static void gunzip(File input, File output) {
        try {
            try (GZIPInputStream in = new GZIPInputStream(new FileInputStream(input));
                 OutputStream out = new FileOutputStream(output)) {

                byte[] buf = new byte[4096];
                int len;
                while ((len = in.read(buf)) > 0) {
                    out.write(buf, 0, len);
                }
            }
        } catch (IOException e){
            throw new GATKException("Exception while unzipping a file:" + input + " to:" + output, e);
        }
    }

    /**
     * Un-gzips the input file to a output file but only if the file's name ends with '.gz'.
     * In this case the new temp file is masked for deletion on exit and returned from this method.
     * Otherwise, that is if the argument is not a gzipped file, this method just returns the argument.
     */
    public static File gunzipToTempIfNeeded(File maybeGzipedFile) {
        if (! maybeGzipedFile.getPath().endsWith(".gz")) {
              return maybeGzipedFile;
        }
        final File result = IOUtils.createTempFile("unzippedFile", "tmp");
        gunzip(maybeGzipedFile, result);
        return result;
    }

    /**
     * Makes a reader for a file, unzipping if the file's name ends with '.gz'.
     */
    public static Reader makeReaderMaybeGzipped(File file) throws IOException {
        final InputStream in = new BufferedInputStream( new FileInputStream(file));
        return makeReaderMaybeGzipped(in, file.getPath().endsWith(".gz"));
    }

    /**
     * makes a reader for an inputStream wrapping it in an appropriate unzipper if necessary
     * @param zipped is this stream zipped
     */
    public static Reader makeReaderMaybeGzipped(InputStream in, boolean zipped) throws IOException {
        if (zipped) {
            return new InputStreamReader(makeZippedInputStream(in));
        } else {
            return new InputStreamReader(in);
        }
    }

    /**
     * creates an input stream from a zipped stream
     * @return tries to create a block gzipped input stream and if it's not block gzipped it produces to a gzipped stream instead
     * @throws ZipException if !in.markSupported()
     */
    public static InputStream makeZippedInputStream(InputStream in) throws IOException {
        Utils.nonNull(in);
        if (BlockCompressedInputStream.isValidFile(in)) {
                return new BlockCompressedInputStream(in);
        } else {
            return new GZIPInputStream(in);
        }
    }

    /**
     * Makes a print stream for a file, gzipping on the fly if the file's name ends with '.gz'.
     */
    public static PrintStream makePrintStreamMaybeGzipped(File file) throws IOException {
        if (file.getPath().endsWith(".gz")) {
            return new PrintStream(new GZIPOutputStream(new FileOutputStream(file)));
        } else {
            return new PrintStream(file);
        }
    }

    /**
     * Creates a temp file that will be deleted on exit
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFile(String name, String extension) {
        try {
            final File file = File.createTempFile(name, extension);
            file.deleteOnExit();

            // Mark corresponding indices for deletion on exit as well just in case an index is created for the temp file:
            new File(file.getAbsolutePath() + Tribble.STANDARD_INDEX_EXTENSION).deleteOnExit();
            new File(file.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION).deleteOnExit();
            new File(file.getAbsolutePath() + ".bai").deleteOnExit();
            new File(file.getAbsolutePath() + ".md5").deleteOnExit();
            new File(file.getAbsolutePath().replaceAll(extension + "$", ".bai")).deleteOnExit();

            return file;
        } catch (IOException ex) {
            throw new GATKException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }


    /**
     * @param extension a file extension, may include 0 or more leading dots which will be replaced with a single dot
     * @return replace the final extension on a path with the given extension
     */
    public static String replaceExtension(String path, String extension){
        Utils.nonNull(path);
        Utils.nonNull(extension);
        final String extensionNoLeadingDot = StringUtils.stripStart(extension, ".");
        return FilenameUtils.removeExtension(path) + '.' + extensionNoLeadingDot;
    }

    public static File replaceExtension(File file, String extension){
        return new File(replaceExtension(file.getPath(), extension));
    }

    /**
     * Schedule a file or directory to be deleted on jvm exit.
     *
     * This will silently delete the directory as well as it's contents.
     * It improves upon {@link FileUtils#forceDeleteOnExit} by deleting directories and files that did not exist at call time.
     * @param dir to be deleted
     */
    public static void deleteRecursivelyOnExit(File dir){
        Runtime.getRuntime().addShutdownHook(new Thread() {

            @Override
            public void run() {
                FileUtils.deleteQuietly(dir);
            }
        });
    }

    /**
     * Converts the given URI to a {@link Path} object. If the filesystem cannot be found in the usual way, then attempt
     * to load the filesystem provider using the thread context classloader. This is needed when the filesystem
     * provider is loaded using a URL classloader (e.g. in spark-submit).
     *
     * @param uriString the URI to convert
     * @return the resulting {@code Path}
     * @throws UserException if an I/O error occurs when creating the file system
     */
    public static Path getPath(String uriString) {
        Utils.nonNull(uriString);
        URI uri = URI.create(uriString);
        try {
            // special case GCS, in case the filesystem provider wasn't installed properly but is available.
            if (CloudStorageFileSystem.URI_SCHEME.equals(uri.getScheme())) {
                return BucketUtils.getPathOnGcs(uriString);
            }
            return uri.getScheme() == null ? Paths.get(uriString) : Paths.get(uri);
        } catch (FileSystemNotFoundException e) {
            try {
                ClassLoader cl = Thread.currentThread().getContextClassLoader();
                if ( cl == null ) {
                    throw e;
                }
                return FileSystems.newFileSystem(uri, new HashMap<>(), cl).provider().getPath(uri);
            }
            catch ( IOException io ) {
                throw new UserException(uriString + " is not a supported path", io);
            }
        }
    }

    /**
     * @param path Path to test
     * @throws org.broadinstitute.hellbender.exceptions.UserException.CouldNotReadInputFile if the file isn't readable
     *         and a regular file
     */
    public static void assertFileIsReadable(final Path path) {
        Utils.nonNull(path);

        try {
            if ( ! Files.exists(path) ) {
                throw new UserException.CouldNotReadInputFile(path, "It doesn't exist.");
            }
            if ( ! Files.isRegularFile(path) ) {
                throw new UserException.CouldNotReadInputFile(path, "It isn't a regular file");
            }
            if ( ! Files.isReadable(path) ) {
                throw new UserException.CouldNotReadInputFile(path, "It is not readable, check the file permissions");
            }
        } catch (com.google.cloud.storage.StorageException cloudBoom) {
            // probably a permissions problem, or perhaps a disabled account.
            // Looks like this for a disabled bucket error:
            //   A USER ERROR has occurred: Couldn't read file gs://foo/bar. Error was:
            //   403: The account for bucket "foo" has been disabled.
            // For no access, it looks like this:
            // (use `gcloud auth application-default revoke` to forget the default credentials)
            //   A USER ERROR has occurred: Couldn't read file gs://(...). Error was:
            //   401: Anonymous users does not have storage.objects.get access to object (...).
            // The user can see the underlying exception by passing
            // -DGATK_STACKTRACE_ON_USER_EXCEPTION=true
            throw new UserException.CouldNotReadInputFile(path, cloudBoom.getCode() + ": " + cloudBoom.getMessage(), cloudBoom);
        }
    }

    /**
     * Checks that one or more user provided files are in fact regular (i.e. not a directory or a special device) readable files.
     *
     * @param files the input files to test.
     * @throws IllegalArgumentException if any input file {@code file} is {@code null} or {@code files} is {@code null}.
     * @throws UserException if any {@code file} is not a regular file or it cannot be read.
     */
    public static void canReadFile( final File... files) {
        Utils.nonNull(files, "Unexpected null input.");
        for (final File file : files) {
            Utils.nonNull(file, "Unexpected null file reference.");
            if (!file.exists()) {
                throw new UserException.CouldNotReadInputFile(file, "The input file does not exist.");
            } else if (!file.isFile()) {
                throw new UserException.CouldNotReadInputFile(file.getAbsolutePath(), "The input file is not a regular file");
            } else if (!file.canRead()) {
                throw new UserException.CouldNotReadInputFile(file.getAbsolutePath(), "The input file cannot be read.  Check the permissions.");
            }
        }
    }

    /**
     * Creates a directory, in local FS, HDFS, or Google buckets to write individual files in.
     */
    public static void createDirectory(final String pathString) throws IOException {
        Utils.nonNull(pathString);

        Files.createDirectory(getPath(pathString));
    }
}
