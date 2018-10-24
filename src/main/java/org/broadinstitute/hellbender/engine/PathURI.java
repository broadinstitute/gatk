package org.broadinstitute.hellbender.engine;

import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
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
