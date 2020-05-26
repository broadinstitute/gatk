package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.FileExtensions;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.file.FileSystems;
import java.nio.file.Path;

/**
 * Interface for htsjdk input/output URIs.
 */
public interface PathURI {

    /**
     * Get a {@code java.net.URI} object for this {@code PathURI}. Will not be null.
     * @return The {@code URI} object for this PathURI.
     */
    URI getURI();

    /**
     * Returns the string from which this {@code PathURI} was originally created. This string may differ from
     * the normalized string returned from a {@code java.nio.file.Path} that has been object resolved from this PathURI.
     * @return string from which this URI as originally created. Will not be null.
     */
    String getURIString();

    /**
     * Return the raw input string that was provided to the constructor.
     */
    String getRawInputString();

    /**
     * @return true if this URI has a scheme that has an installed {@code java.nio} file system provider
     * ({@linktourl https://docs.oracle.com/javase/8/docs/api/java/nio/file/spi/FileSystemProvider.html}). This does not
     * guarantee the URI can be converted into a {@code java.nio.file.Path}, since the URI can be syntactically
     * valid, and specify a valid file system provider, but still fail to be semantically meaningful.
     */
    boolean hasFileSystemProvider();

    /**
     * Return true if this {code PathURI} can be resolved to an {@code java.nio} Path. If true, {@code #toPath()} can be
     * safely called.
     *
     * There are cases where a valid URI with a valid scheme backed by an installed {@code java.nio File System
     * still can't be turned into a {@code java.nio.file.Path}, i.e., the following specifies an invalid
     * authority "namenode":
     *
     *  hdfs://namenode/to/file
     *
     * The current implementation returns false for these cases (toPath will fail, getInvalidPathReason
     * returns the reason code).
     */
    boolean isPath();

    /**
     * Resolve this PathURI to an NIO {@code java.nio.file.Path}. Can be safely called only if {@code #isPath()} returns true.
     */
    Path toPath();

    /**
     * Return a string message describing why this URI cannot be converted to a {@code java.nio.file.Path}
     * ({@code #isPath()} returns false).
     * @return Message explaining toPath failure reason, since it can fail for various reasons.
     */
    String getToPathFailureReason();

    /**
     * Return the scheme for this PathURI. For file URIs (URIs that have no explicit scheme), this will return
     * the scheme "file".
     * @return the scheme String or this URI, if any. May be null.
     */
    default String getScheme() {
        return getURI().getScheme();
    }

    /**
     * @return the extension of the last component of the hierarchical part of the scheme-specific part of the
     * URI, if any, including the ".". Note that this only returns the part of the last component after the last
     * ".", ie. it will return ".gz" for a name that ends in ".fasta.gz" (the {@link #hasExtension(String)} method
     * can be used to test for the presence of multi-part extensions such as this). If the hierarchical name ends
     * with a last component that does not contain a ".", returns an empty String.
     * @throws IllegalArgumentException if the hierarchical name ends with the default file system separator
     * (i.e. "/").
     */
    default String getExtension() {
        final String hierarchicalPath = getURI().getPath();
        final int indexOfLastComponent = hierarchicalPath.lastIndexOf(FileSystems.getDefault().getSeparator());
        if (indexOfLastComponent != -1 && indexOfLastComponent < hierarchicalPath.length() - 1) {
            final String lastComponent = hierarchicalPath.substring(indexOfLastComponent + 1);
            if (lastComponent.length() > 0) {
                final int indexOfLastDot = lastComponent.lastIndexOf('.');
                if (indexOfLastDot != -1 && indexOfLastDot < lastComponent.length() - 1) {
                    // return a string that includes the leading "." to enable easy comparison with the many
                    // internal file extension constants we have that include the leading "." (i.e., in htsjdk),
                    // and also for API consistency (since hasExtension() requires the candidate extension to
                    // include a leading ".", this allows hasExtension(getExtension()) to always work whenever
                    // getExtension() succeeds)
                    return lastComponent.substring(indexOfLastDot);
                }
            }
            return "";
        }
        throw new IllegalArgumentException(
                String.format(
                        "Input (%s) contains only a path, with no name component from which to retrieve an extension", this));
    }

    /**
     * Return true if the path component (the hierarchical part of the scheme specific part of the underlying URI)
     * ends with the provided {@code extension} string. This method can be used to test for both single and
     * multi-part extensions (ie. for a name that ends in ".fasta.gz", it will return true for both ".gz" and
     * ".fasta.gz".
     *
     * @param extension the target extension to test, INCLUDING the leading ".". May not be null.
     * @return true if the path component of this specifier ends with the extension, otherwise false.
     */
    default boolean hasExtension(final String extension) {
        Utils.nonNull(extension, "Target extension must not be null");
        Utils.validateArg(extension.length() > 0, "Target extension must be length > 0");
        Utils.validateArg(extension.charAt(0) == '.', "Target extension must include the leading '.'");

        // We don't want to use {@code #getExtension} here, since it won't work correctly if we're comparing an
        // extension that uses multiple . chars, such as .fasta.gz.
        return getURI().getPath().toLowerCase().endsWith(extension.toLowerCase());
    }

    /**
     * @return the base name (the last component of the hierarchical part of the scheme-specific part of the URI,
     * after the last "/"), up to but not including the extension (the last ".").
     * @throws IllegalArgumentException if the last component is empty (ie, the component ends in "/"), or the last
     * component exists but starts with "."
     */
    default String getBaseName() {
        final String hierarchicalPath = getURI().getPath();
        final int indexOfLastComponent = hierarchicalPath.lastIndexOf(FileSystems.getDefault().getSeparator());
        if (indexOfLastComponent != -1 && indexOfLastComponent < hierarchicalPath.length() - 1) {
            final String lastComponent = hierarchicalPath.substring(indexOfLastComponent + 1);
            if (lastComponent.length() > 0) {
                final int indexOfLastDot = lastComponent.lastIndexOf('.');
                if (indexOfLastDot == -1) {
                    return lastComponent;
                } else if (indexOfLastDot > 1) {
                    return lastComponent.substring(0, indexOfLastDot);
                }
            }
            return "";
        }
        throw new IllegalArgumentException(String.format("Input (%s) contains only a path, but no name component", this));
    }

    /**
     * Returns true if the file's extension is ".sam"".
     */
    default boolean isSam() {
        return hasExtension(FileExtensions.SAM);
    }

    /**
     * Returns true if the file's extension is ".bam"".
     */
    default boolean isBam() {
        return hasExtension(FileExtensions.BAM);
    }

    /**
     * Returns true if the GATKPathSpecifier's extension is ".cram".
     */
    default boolean isCram() {
        return hasExtension(FileExtensions.CRAM);
    }

    /**
     * @return true if this path spec has a FASTA file extension
     */
    default boolean isFasta() {
        return FileExtensions.FASTA.stream().anyMatch(this::hasExtension);
    }

    /**
     * Get a {@code InputStream} for this URI if a provider is for the URI's scheme is available.
     * @return {@code InputStream} for this URI.
     */
    InputStream getInputStream();

    /**
     * Get an {@code OutputStream} for this URI if a provider is for the URI's scheme.
     * @return {@code OutputStream} for this URI.
     */
    OutputStream getOutputStream();
}
