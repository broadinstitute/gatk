package org.broadinstitute.hellbender.engine;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;

import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * ArgumentsBuilder specialization that includes management of a temporary JIMFS file system. This class
 * can be used to test Path URIs on inputs that live on a java.nio file system, and are not java.io.File
 * compatible.
 *
 * Methods accept File arguments, which are copied to a Jimfs file system managed by instacnes o this class,
 * and Path URIs specifying the Jimfs copy are used as the actual argument values.
 *
 * TODO: Test methods need to declare IOExceptions since many of these methods throw IOException directly?
 * TODO: This has to live in the test hierarchy in order to be able to import jimfs. Not sure why ArgumentsBuilder isn't also here
 *
 * Usage is intended to be via try-with-resources:
 *
 *     try (final ArgumentsBuilderWithJIMFS args = new ArgumentsBuilderWithJIMFS()) {
 *         // add an input file with corresponding index
 *         args.addFileWithIndexFromJIMFS(StandardArgumentDefinitions.INPUT_SHORT_NAME, cramFile, indexFile);
 *
 *         // add the reference, copying the reference, companion dictionary, and index to jimfs
 *         args.addReferenceFileFromJIMFS(referenceFile);
 *
 *         final File outputFile = createTempFile("testReadCramWithIndexFromPath", ".bam");
 *         args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
 *         args.add(outputFile.getAbsolutePath());
 *
 *         runCommandLine(args);
 *
 *         // Assert.assert...
 *     }
 */
public class ArgumentsBuilderWithJIMFS extends ArgumentsBuilder implements Closeable {
    final FileSystem jimfs;
    final Path jimfsRootPath;
    final Path tempDir;
    final static String TEMP_DIR = "/argsBuilderJIMFSTemp";

    // Create a new in-memory jimfs file system
    public ArgumentsBuilderWithJIMFS() throws IOException {
        jimfs = Jimfs.newFileSystem(Configuration.unix());
        jimfsRootPath = jimfs.getRootDirectories().iterator().next();
        tempDir = jimfsRootPath.resolve(TEMP_DIR);
        Files.createDirectory(tempDir);
    }

    @Override
    public void close() throws IOException {
        jimfs.close();
    }

    /**
     * add a File as a Path/URI by first copying it to JIMFS and then using the resulting Path as the arg value
     */
    public ArgumentsBuilder addInputFromJIMFS(final File inputFile) throws IOException {
        Utils.nonNull(inputFile);

        add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        final Path path = copyAndResolve(inputFile);
        add(path.toUri().toString());

        return this;
    }

    /**
     * add a File as a Path/URI by first copying it and it's companion index file to JIMFS and then using the
     * resulting the resulting Path URI as the argument value
     */
    public ArgumentsBuilder addInputWithIndexFromJIMFS(final File inputFile, final File indexFile) throws IOException {
        Utils.nonNull(inputFile);
        Utils.nonNull(indexFile);

        add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        final Path path = copyAndResolve(inputFile);
        copyAndResolve(indexFile);
        add(path.toUri().toString());

        return this;
    }

    /**
     * add a File as a Path/URI by first copying it to JIMFS and then using the resulting Path URI
     * as the argument value
     */
    public ArgumentsBuilder addReferenceFromJIMFS(final File referenceFile) throws IOException {
        Utils.nonNull(referenceFile);

        // copy the reference, reference dictionary, and reference index
        final Path path = copyAndResolve(referenceFile);
        copyAndResolve(ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(referenceFile));
        copyAndResolveToName(ReferenceSequenceFileFactory.getFastaIndexFileName(referenceFile.toPath()).toFile(),
                referenceFile.getName() + ".fai");

        add("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        add(path.toUri().toString());

        return this;
    }

    // TODO: This violates the builder pattern since it needs to return the Path so the test can use it to validate output.
    public Path addOutputFromJIMFS(final String requestedName) throws IOException {
        final Path tempOut = getTempFileFromJIMFS(requestedName);
        add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        add(tempOut.toUri().toString());
        return tempOut;
    }

    /**
     * Obtain a temporary Path to be used as an output file.
     * @param requestedName name of the requested output file
     * @return Path of the newly created file in the jimfs temp directory
     */
    public Path getTempFileFromJIMFS(final String requestedName) throws IOException {
        Path tempFile = tempDir.resolve(requestedName);
        if (Files.exists(tempFile)) {
            throw new IllegalArgumentException("The requested temporary file name has already been used");
        } else {
            return IOUtils.getPath(tempFile.toUri().toString());
        }
    }

    private Path copyAndResolve(final File inputFile) throws IOException {
        Utils.nonNull(inputFile);
        return copyAndResolveToName(inputFile, inputFile.getName());
    }

    private Path copyAndResolveToName(final File inputFile, final String targetName) throws IOException {
        Utils.nonNull(inputFile);
        return Files.copy(inputFile.toPath(), jimfsRootPath.resolve(targetName));
    }
}
