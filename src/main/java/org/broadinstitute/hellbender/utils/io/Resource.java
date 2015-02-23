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

import java.io.File;
import java.io.InputStream;

/**
 * Stores a resource by path and a relative class.
 */
public class Resource {
    private final String path;
    private final Class<?> relativeClass;
    private final ClassLoader relativeClassLoader;

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
        ClassLoader classLoader = null;
        if (relativeClass != null)
            classLoader = relativeClass.getClassLoader();
        this.relativeClassLoader = classLoader != null ? classLoader : ClassLoader.getSystemClassLoader();
    }

    public Class<?> getRelativeClass() {
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
}
