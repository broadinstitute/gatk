package org.broadinstitute.hellbender.utils.io;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.util.TabixUtils;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.GetSampleName;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;

import java.io.*;
import java.net.URI;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

public final class IOUtils {
    private static final Logger logger = LogManager.getLogger(IOUtils.class);
    private static final File DEV_DIR = new File("/dev");

    // see https://support.hdfgroup.org/HDF5/doc/H5.format.html
    private static final byte hdf5HeaderSignature[] = { (byte) 0x89, 'H', 'D', 'F', '\r', '\n', (byte) 0x1A, '\n' };

    /**
     * Returns true if the file's extension is CRAM.
     */
    public static boolean isCramFile(final File inputFile) {
        return isCramFileName(inputFile.getName());
    }

    /**
     * Returns true if the file's extension is CRAM.
     */
    public static boolean isCramFile(final Path path) {
        return isCramFileName(path.getFileName().toString());
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
     * Given a Path, determine if it is an HDF5 file without requiring that we're on a platform that supports
     * HDF5 (let the caller decide if a return value of false is fatal).
     *
     * @param hdf5Candidate a Path representing the input to be inspected
     * @return true if the candidate Path is an HDF5 file, otherwise false
     */
    public static boolean isHDF5File(final Path hdf5Candidate) {
        try (final DataInputStream candidateStream = new DataInputStream(Files.newInputStream(hdf5Candidate))) {
            final byte candidateHeader[] = new byte[hdf5HeaderSignature.length];
            candidateStream.read(candidateHeader, 0, candidateHeader.length);
            return Arrays.equals(candidateHeader, hdf5HeaderSignature);
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(String.format("I/O error reading from input stream %s", hdf5Candidate), e);
        }
    }

    /**
     * Creates a temp directory with the given prefix.
     *
     * The directory and any contents will be automatically deleted at shutdown.
     *
     * This will not work if the temp dir is not representable as a File.
     *
     * @param prefix       Prefix for the directory name.
     * @return The created temporary directory.
     */
    public static File createTempDir(String prefix) {
        try {
            final File tmpDir = Files.createTempDirectory(prefix)
                    .normalize()
                    .toFile();
            deleteRecursivelyOnExit(tmpDir);

            return tmpDir;
        } catch (final IOException | SecurityException e) {
            throw new UserException.BadTempDir(e.getMessage(), e);
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
            File tempFile = File.createTempFile(prefix, suffix, directory).toPath().normalize().toFile();
            FileUtils.writeStringToFile(tempFile, content, StandardCharsets.UTF_8);
            return tempFile;
        } catch (IOException e) {
            throw new UserException.BadTempDir(e.getMessage(), e);
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
     * Writes an embedded resource to a temporary file. The temporary file is automatically scheduled for deletion
     * on exit.
     * @param resource Embedded resource.
     * @return the temporary file containing the contents of the resource, which is automatically scheduled for
     * deletion on exit.
     */
    public static File writeTempResource(final Resource resource) {
        final File tempFile = createTempFile(
                FilenameUtils.getBaseName(resource.getPath()) + ".",
                "." + FilenameUtils.getExtension(resource.getPath()));
        writeResource(resource, tempFile);
        return tempFile;
    }

    /**
     * Create a resource from a path and a relative class, and write it to a temporary file.
     * If the relative class is null then the system classloader will be used and the path must be absolute.
     * The temporary file is automatically scheduled for deletion on exit.
     * @param resourcePath Relative or absolute path to the class.
     * @param relativeClass Relative class to use as a class loader and for a relative package.
     * @return a temporary file containing the contents of the resource, which is automatically scheduled
     * for deletion on exit.
     */
    public static File writeTempResourceFromPath(final String resourcePath, final Class<?> relativeClass) {
        Utils.nonNull(resourcePath, "A resource path must be provided");
        final Resource resource = new Resource(resourcePath, relativeClass);
        return writeTempResource(resource);
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
     * Extracts the tar.gz file given by {@code tarGzFilePath}.
     * Input {@link Path} MUST be to a gzipped tar file.
     * Will extract contents in the containing folder of {@code tarGzFilePath}.
     * Will throw an exception if files exist already.
     * @param tarGzFilePath {@link Path} to a gzipped tar file for extraction.
     */
    public static void extractTarGz(final Path tarGzFilePath) {
        extractTarGz(tarGzFilePath, tarGzFilePath.getParent(), false);
    }

    /**
     * Extracts the tar.gz file given by {@code tarGzFilePath}.
     * Input {@link Path} MUST be to a gzipped tar file.
     * Will throw an exception if files exist already.
     * @param tarGzFilePath {@link Path} to a gzipped tar file for extraction.
     * @param destDir {@link Path} to the directory where the contents of {@code tarGzFilePath} will be extracted.
     */
    public static void extractTarGz(final Path tarGzFilePath, final Path destDir) {
        extractTarGz(tarGzFilePath, destDir, false);
    }

    /**
     * Extracts the tar.gz file given by {@code tarGzFilePath}.
     * Input {@link Path} MUST be to a gzipped tar file.
     * @param tarGzFilePath {@link Path} to a gzipped tar file for extraction.
     * @param destDir {@link Path} to the directory where the contents of {@code tarGzFilePath} will be extracted.
     * @param overwriteExistingFiles If {@code true}, will enable overwriting of existing files.  If {@code false}, will cause an exception to be thrown if files exist already.
     */
    public static void extractTarGz(final Path tarGzFilePath, final Path destDir, final boolean overwriteExistingFiles) {

        logger.info("Extracting data from archive: " + tarGzFilePath.toUri());

        // Create a stream for the data sources input.
        // (We know it will be a tar.gz):
        try ( final InputStream fi = Files.newInputStream(tarGzFilePath);
              final InputStream bi = new BufferedInputStream(fi);
              final InputStream gzi = new GzipCompressorInputStream(bi);
              final TarArchiveInputStream archiveStream = new TarArchiveInputStream(gzi)) {

            extractFilesFromArchiveStream(archiveStream, tarGzFilePath, destDir, overwriteExistingFiles);
        }
        catch (final IOException ex) {
            throw new UserException("Could not extract data from: " + tarGzFilePath.toUri(), ex);
        }
    }

    private static void extractFilesFromArchiveStream(final TarArchiveInputStream archiveStream,
                                                      final Path localTarGzPath,
                                                      final Path destDir,
                                                      final boolean overwriteExistingFiles) throws IOException {

        // Adapted from: http://commons.apache.org/proper/commons-compress/examples.html

        // Go through the archive and get the entries:
        TarArchiveEntry entry;
        while ((entry = archiveStream.getNextTarEntry()) != null) {

            logger.info("Extracting file: " + entry.getName());

            // Make sure we can read the data for the entry:
            if (!archiveStream.canReadEntryData(entry)) {
                throw new UserException("Could not read data from archive file(" + localTarGzPath.toUri() + "): " + entry.getName());
            }

            // Get the path for the entry on disk and make sure it's OK:
            final Path extractedEntryPath = destDir.resolve(entry.getName()).normalize();
            ensurePathIsOkForOutput(extractedEntryPath, overwriteExistingFiles);

            // Now we can create the entry in our output location:
            if (entry.isDirectory()) {
                Files.createDirectories(extractedEntryPath);
            }
            else {
                // Make sure the parent directory exists:
                Files.createDirectories(extractedEntryPath.getParent());

                if ( entry.isFIFO() ) {
                    // Handle a fifo file:
                    createFifoFile(extractedEntryPath, overwriteExistingFiles);
                }
                else if ( entry.isSymbolicLink() ) {
                    // Handle a symbolic link:
                    final String linkName = entry.getLinkName();

                    // If the link already exists, we must clear it:
                    if ( Files.exists(extractedEntryPath) && overwriteExistingFiles ) {
                        removeFileWithWarning(extractedEntryPath);
                    }

                    Files.createSymbolicLink(extractedEntryPath, Paths.get(linkName));
                }
                else if ( entry.isLink() ) {
                    // Handle a hard link:
                    final String linkName = entry.getLinkName();

                    // If the link already exists, we must clear it:
                    if ( Files.exists(extractedEntryPath) && overwriteExistingFiles ) {
                        removeFileWithWarning(extractedEntryPath);
                    }

                    Files.createLink(extractedEntryPath, Paths.get(linkName));
                }
                else if ( entry.isFile() ) {
                    // Handle a (default) file entry:

                    // Create the output file from the stream:
                    try (final OutputStream o = Files.newOutputStream(extractedEntryPath)) {
                        org.apache.commons.io.IOUtils.copy(archiveStream, o);
                    }
                }
                else {
                    // Right now we don't know how to handle any other file types:
                    throw new UserException("Cannot extract file from tar.gz (unknown type): " + entry.toString());
                }
            }
        }
    }

    private static void ensurePathIsOkForOutput(final Path p, final boolean overwriteExistingFiles) {
        if ( Files.exists(p) ) {
            if ( overwriteExistingFiles ) {
                logger.warn("Overwriting existing output destination: " + p.toUri());
            }
            else {
                throw new UserException("Output destination already exists: " + p.toUri());
            }
        }
    }

    /**
     * Create a Unix FIFO file with the given path string.
     * If requested file already exists, will throw an exception.
     * Will throw an Exception on failure.
     * @param fifoFilePath {@link Path} to the FIFO file to be created.
     * @return The {@link File} object pointing to the created FIFO file.
     */
    public static File createFifoFile(final Path fifoFilePath) {
        return createFifoFile(fifoFilePath, false);
    }

    private static void removeFileWithWarning(final Path filePath) {
        logger.warn("File already exists in path.  Replacing existing file: " + filePath.toUri());
        try {
            Files.delete(filePath);
        }
        catch (final IOException ex) {
            throw new UserException("Could not replace existing file: " + filePath.toUri());
        }
    }

    /**
     * Create a Unix FIFO file with the given path string.
     * Will throw an Exception on failure.
     * @param fifoFilePath {@link Path} to the FIFO file to be created.
     * @param overwriteExisting If {@code true} will overwrite an existing file in the requested location for the FIFO file.  If {@code false} will throw an exception if the file exists.
     * @return The {@link File} object pointing to the created FIFO file.
     */
    public static File createFifoFile(final Path fifoFilePath, final boolean overwriteExisting) {

        // Make sure we're allowed to create the file:
        if ( Files.exists(fifoFilePath) ) {
            if ( (!overwriteExisting) ) {
                throw new UserException("Cannot create fifo file.  File already exists: " + fifoFilePath.toUri());
            }
            else {
                removeFileWithWarning(fifoFilePath);
            }
        }

        // Create the FIFO by executing mkfifo via another ProcessController
        final ProcessSettings mkFIFOSettings = new ProcessSettings(new String[]{"mkfifo", fifoFilePath.toFile().getAbsolutePath()});
        mkFIFOSettings.getStdoutSettings().setBufferSize(-1);
        mkFIFOSettings.setRedirectErrorStream(true);

        // Now perform the system call:
        final ProcessController mkFIFOController = new ProcessController();
        final ProcessOutput     result           = mkFIFOController.exec(mkFIFOSettings);
        final int               exitValue        = result.getExitValue();

        final File fifoFile = fifoFilePath.toFile();

        // Make sure we're OK:
        if (exitValue != 0) {
            throw new GATKException(String.format(
                    "Failure creating FIFO named (%s). Got exit code (%d) stderr (%s) and stdout (%s)",
                    fifoFilePath.toFile().getAbsolutePath(),
                    exitValue,
                    result.getStderr() == null ? "" : result.getStderr().getBufferString(),
                    result.getStdout() == null ? "" : result.getStdout().getBufferString()));
        } else if (!fifoFile.exists()) {
            throw new GATKException(String.format("FIFO (%s) created but doesn't exist", fifoFilePath.toFile().getAbsolutePath()));
        } else if (!fifoFile.canWrite()) {
            throw new GATKException(String.format("FIFO (%s) created isn't writable", fifoFilePath.toFile().getAbsolutePath()));
        }

        return fifoFile;
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
        return createTempFileInDirectory(name, extension, null);
    }

    /**
     * Creates a temp file in a target directory that will be deleted on exit
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file name.
     * @param targetDir Directory in which to create the temp file. If null, the default temp directory is used.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFileInDirectory(final String name, String extension, final File targetDir) {
        try {

            if ( !extension.startsWith(".") ) {
                extension = "." + extension;
            }

            final File file = File.createTempFile(name, extension, targetDir);
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
     * Creates a temp path that will be deleted on exit.
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     *
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     *
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static Path createTempPath(String name, String extension) {
        try {

            if ( !extension.startsWith(".") ) {
                extension = "." + extension;
            }

            final Path path = Files.createTempFile(getPath(System.getProperty("java.io.tmpdir")), name, extension);
            IOUtil.deleteOnExit(path);

            // Mark corresponding indices for deletion on exit as well just in case an index is created for the temp file:
            final String filename = path.getFileName().toString();
            IOUtil.deleteOnExit(path.resolveSibling(filename + Tribble.STANDARD_INDEX_EXTENSION));
            IOUtil.deleteOnExit(path.resolveSibling(filename + TabixUtils.STANDARD_INDEX_EXTENSION));
            IOUtil.deleteOnExit(path.resolveSibling(filename + BAMIndex.BAMIndexSuffix));
            IOUtil.deleteOnExit(path.resolveSibling(filename.replaceAll(extension + "$", ".bai")));
            IOUtil.deleteOnExit(path.resolveSibling(filename + ".md5"));

            return path;
        } catch (final IOException ex) {
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
     * Also makes an attempt to interpret the argument as a file name if it's not a URI.
     *
     * @param uriString the URI to convert.
     * @return the resulting {@code Path}
     * @throws UserException if an I/O error occurs when creating the file system
     */
    public static Path getPath(String uriString) {
        Utils.nonNull(uriString);
        URI uri;
        try {
            uri = URI.create(uriString);
        } catch (IllegalArgumentException x) {
            // not a valid URI. Caller probably just gave us a file name.
            return Paths.get(uriString);
        }
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
            catch (ProviderNotFoundException x) {
                // not a valid URI. Caller probably just gave us a file name or "chr1:1-2".
                return Paths.get(uriString);
            }
            catch ( IOException io ) {
                throw new UserException(uriString + " is not a supported path", io);
            }
        }
    }

    /**
     * Gets the absolute Path name with the URI marker, handling the special case of the default file system by removing
     * the file:// prefix.
     *
     * @param path path to get the absolute name.
     * @return a String with the absolute name, and the file:// protocol removed, if it was present.
     */
    public static String getAbsolutePathWithoutFileProtocol(final Path path) {
        return path.toAbsolutePath().toUri().toString().replaceFirst("^file://", "");
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

    public static final String urlEncode(final String string) {
        try {
            return URLEncoder.encode(string, GetSampleName.STANDARD_ENCODING);
        } catch (final UnsupportedEncodingException ex) {
            throw new UserException("Could not encode sample name", ex);
        }
    }

    public static final String urlDecode(final String string) {
        try {
            return URLDecoder.decode(string, GetSampleName.STANDARD_ENCODING);
        } catch (final UnsupportedEncodingException ex) {
            throw new UserException("Could not decode sample name", ex);
        }
    }
}
