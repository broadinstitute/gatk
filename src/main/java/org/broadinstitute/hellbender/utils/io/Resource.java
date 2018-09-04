package org.broadinstitute.hellbender.utils.io;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

/**
 * Stores a resource by path and a relative class.
 */
public final class Resource {
    private final String path;
    private final Class<?> relativeClass;

    /**
     * Path to the directory for large runtime system resources
     */
    public final static String LARGE_RUNTIME_RESOURCES_PATH = "large";

    /**
     * Create a resource with a path and a relative class.
     * @param path Relative or absolute path to the class.
     * @param relativeClass Relative class to use as a class loader and for a relative package.
     *
     * If the relative class is null then the system classloader will be used and the path must be absolute.
     */
    public Resource(String path, Class<?> relativeClass) {
        this.path = path;
        this.relativeClass = relativeClass;
    }

    public final Class<?> getRelativeClass() {
        return relativeClass;
    }

    public String getPath() {
        return path;
    }

    public String getFullPath() {
        if (relativeClass == null)
            return path;
        if (new File(path).isAbsolute())
            return path;
        return String.format("%s%s%s",
                relativeClass.getPackage().getName().replace('.', File.separatorChar),
                File.separator,
                path);
    }

    /**
     * Get the contents of this resource as an InputStream
     * @throws IllegalArgumentException if resource cannot be read
     * @return an input stream that will read the contents of this resource
     */
    public InputStream getResourceContentsAsStream() {
        final Class<?> clazz = getRelativeClass();

        final InputStream inputStream;
        if (clazz == null) {
            inputStream = ClassLoader.getSystemResourceAsStream(path);
            if (inputStream == null)
                throw new IllegalArgumentException("Resource not found: " + path);
        } else {
            inputStream = clazz.getResourceAsStream(path);
            if (inputStream == null)
                throw new IllegalArgumentException("Resource not found relative to " + clazz + ": " + path);

        }

        return inputStream;
    }

    /**
     *  Given a resource path that is either on the file system or in a jar, return a File.
     *
     *  Note:  A side-effect is that the temp file will be copied to temporary directory.
     *
     * @param resourcePath location of the resource.  If in a jar file, use relative path.
     * @return the File in a temp directory
     * @throws IOException
     */
    public static File getResourceContentsAsFile(final String resourcePath) throws IOException {
        final File tmpResourceFile = File.createTempFile("tmp_read_resource_", ".config");
        final InputStream systemResourceAsStream = ClassLoader.getSystemResourceAsStream(resourcePath);

        if (systemResourceAsStream == null) {
            throw new GATKException("Null value when trying to read system resource.  Cannot find: " + resourcePath);
        }

        FileUtils.copyInputStreamToFile(systemResourceAsStream,
                tmpResourceFile);
        return tmpResourceFile;
    }
}
