package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.*;
import java.nio.file.spi.FileSystemProvider;

/**
 * Default implementation for PathURI.
 *
 * This class takes a raw string that is to be interpreted as a path specifier, and converts it internally to a
 * URI and/or Path object. If no scheme is provided as part of the raw string used in the constructor(s), the
 * input is assumed to represent a file on the local file system, and will be backed by a URI with a "file:/"
 * scheme and a path part that is automatically encoded/escaped to ensure it is a valid URI. If the raw string
 * contains a scheme, it will be backed by a URI formed from the raw string as presented, with no further
 * encoding/escaping.
 *
 * For example, a URI that contains a scheme and has an embedded "#" in the path will be treated as a URI
 * having a fragment delimiter. If the URI contains an scheme, the "#" will be escaped and the encoded "#"
 * will be interpreted as part of the URI path.
 *
 * There are 3 succeeding levels of input validation/conversion:
 *
 * 1) PathSpecifier constructor: requires a syntactically valid URI, possibly containing a scheme (if no scheme
 *    is present the path part will be escaped/encoded), or a valid local file reference
 * 2) hasFileSystemProvider: true if the input string is an identifier that is syntactically valid, and is backed by
 *    an installed {@code java.nio} file system provider that matches the URI scheme
 * 3) isPath: syntactically valid URI that can be resolved to a java.io.Path by the associated provider
 *
 * Definitions taken from RFC 2396 "Uniform Resource Identifiers (URI): Generic Syntax"
 * (https://www.ietf.org/rfc/rfc2396.txt):
 *
 * "absolute" URI  - specifies a scheme
 * "relative" URI  - does not specify a scheme
 * "opaque" URI - an "absolute" URI whose scheme-specific part does not begin with a slash character
 * "hierarchical" URI - either an "absolute" URI whose scheme-specific part begins with a slash character,
 *  or a "relative" URI (no scheme)
 *
 * URIs that do not make use of the slash "/" character for separating hierarchical components are
 * considered "opaque" by the generic URI parser.
 *
 * General syntax for an "absolute" URI:
 *
 *     <scheme>:<scheme-specific-part>
 *
 * Many "hierarchical" URI schemes use this syntax:
 *
 *     <scheme>://<authority><path>?<query>
 *
 * More specifically:
 *
 *     absoluteURI   = scheme ":" ( hier_part | opaque_part )
 *         hier_part     = ( net_path | abs_path ) [ "?" query ]
 *         net_path      = "//" authority [ abs_path ]
 *         abs_path      = "/"  path_segments
 *         opaque_part   = uric_no_slash *uric
 *         uric_no_slash = unreserved | escaped | ";" | "?" | ":" | "@" | "&" | "=" | "+" | "$" | ","
 */
public class PathSpecifier implements PathURI, Serializable {
    private static final long serialVersionUID = 1L;

    private final String    rawInputString;     // raw input string provided by th user; may or may not have a scheme
    private final URI       uri;                // working URI; always has a scheme (assume "file" if not provided)
    private transient String pathFailureReason; // cache the reason for "toPath" conversion failure
    private transient Path  cachedPath;         // cache the Path associated with this URI if its "Path-able"

    /**
     * If the raw input string already contains a scheme (including a "file" scheme), assume its already
     * properly escape/encoded. If no scheme component is present, assume it referencess a raw path on the
     * local file system, so try to get a Path first, and then retrieve the URI from the resulting Path.
     * This ensures that input strings that are local file references without a scheme component and contain
     * embedded characters are valid in file names, but which would otherwise be interpreted as excluded
     * URI characters (such as the URI fragment delimiter "#") are properly escape/encoded.
     * @param rawInputString
     */
    public PathSpecifier(final String rawInputString) {
        Utils.nonNull(rawInputString);
        this.rawInputString = rawInputString;

        URI tempURI;
        try {
            tempURI = new URI(rawInputString);
            if (!tempURI.isAbsolute()) {
                // if the URI has no scheme, assume its a local (non-URI) file reference, and resolve
                // it to a Path and retrieve the URI from the Path to ensure proper escape/encoding
                setCachedPath(Paths.get(rawInputString));
                tempURI = getCachedPath().toUri();
            }
        } catch (URISyntaxException uriException) {
            // the input string isn't a valid URI; assume its a local (non-URI) file reference, and
            // use the URI resulting from the corresponding Path
            try {
                setCachedPath(Paths.get(rawInputString));
                tempURI = getCachedPath().toUri();
            } catch (InvalidPathException | UnsupportedOperationException | SecurityException pathException) {
                // we have two exceptions, each of which might be relevant since we can't tell whether
                // the user intended to provide a local file reference or a URI, so preserve both
                final String errorMessage = String.format(
                        "%s can't be interpreted as a local file (%s) or as a URI (%s).",
                        rawInputString,
                        pathException.getMessage(),
                        uriException.getMessage());
                throw new IllegalArgumentException(errorMessage, pathException);
            }
        }
        if (!tempURI.isAbsolute()) {
            // assert the invariant that every URI we create has a scheme, even if the raw input string does not
            throw new UserException("URI has no scheme");
        }

        uri = tempURI;
    }

    @Override
    public URI getURI() {
        return uri;
    }

    @Override
    public String getURIString() {
        return getURI().toString();
    }

    /**
     * Return the raw input string provided to the constructor.
     */
    @Override
    public String getRawInputString() { return rawInputString; }

    @Override
    public boolean hasFileSystemProvider() {
        // try to find a provider; assume that our URI always has a scheme
        for (FileSystemProvider provider: FileSystemProvider.installedProviders()) {
            if (provider.getScheme().equalsIgnoreCase(uri.getScheme())) {
                return true;
            }
        }
        return false;
    }

    @Override
    public boolean isPath() {
        try {
            return getCachedPath() != null || toPath() != null;
        } catch (ProviderNotFoundException |
                FileSystemNotFoundException |
                IllegalArgumentException |
                UserException |
                AssertionError e) {
            // jimfs throws an AssertionError that wraps a URISyntaxException when trying to create path where
            // the scheme-specific part is missing or incorrect
            pathFailureReason = e.getMessage();
            return false;
        }
    }

    /**
     * Resolve the URI to a {@link Path} object.
     *
     * @return the resulting {@code Path}
     * @throws UserException if an I/O error occurs when creating the file system
     */
    @Override
    public Path toPath() {
        if (getCachedPath() != null) {
            return getCachedPath();
        } else {
            final Path tmpPath = Paths.get(getURI());
            setCachedPath(tmpPath);
            return tmpPath;
        }
    }

    @Override
    public String getToPathFailureReason() {
        if (pathFailureReason == null) {
            try {
                toPath();
                return String.format("'%s' appears to be a valid Path", rawInputString);
            } catch (ProviderNotFoundException e) {
                return String.format("ProviderNotFoundException: %s", e.getMessage());
            } catch (FileSystemNotFoundException e) {
                return String.format("FileSystemNotFoundException: %s", e.getMessage());
            } catch (IllegalArgumentException e) {
                return String.format("IllegalArgumentException: %s", e.getMessage());
            } catch (UserException e) {
                return String.format("UserException: %s", e.getMessage());
            }
        }
        return pathFailureReason;
    }

    @Override
    public InputStream getInputStream() {
        if (!isPath()) {
            throw new UserException(getToPathFailureReason());
        }

        final Path resourcePath = toPath();
        try {
            return Files.newInputStream(resourcePath);
        } catch (IOException e) {
            throw new UserException(
                    String.format("Could not create open input stream for %s (as URI %s)", getRawInputString(), getURIString()), e);
        }
    }

    @Override
    public OutputStream getOutputStream() {
        if (!isPath()) {
            throw new UserException(getToPathFailureReason());
        }

        final Path resourcePath = toPath();
        try {
            return Files.newOutputStream(resourcePath);
        } catch (IOException e) {
            throw new UserException(String.format("Could not open output stream for %s (as URI %s)", getRawInputString(), getURIString()), e);
        }
    }

    // get the cached path associated with this URI if its already been created
    protected Path getCachedPath() { return cachedPath; }

    protected void setCachedPath(Path path) {
        this.cachedPath = path;
    }

    @Override
    public String toString() {
        return rawInputString;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PathSpecifier)) return false;

        PathSpecifier that = (PathSpecifier) o;

        if (!getRawInputString().equals(that.getRawInputString())) return false;
        if (!getURI().equals(that.getURI())) return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = getRawInputString().hashCode();
        result = 31 * result + getURI().hashCode();
        return result;
    }

}

