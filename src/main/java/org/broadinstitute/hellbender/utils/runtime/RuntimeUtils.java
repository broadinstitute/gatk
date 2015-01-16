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

package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.lang3.StringUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class RuntimeUtils {
    public static final String[] PATHS;

    static {
        String path = System.getenv("PATH");
        if (path == null)
            path = System.getenv("path");
        if (path == null) {
            PATHS = new String[0];
        } else {
            PATHS = StringUtils.split(path, File.pathSeparatorChar);
        }
    }

    /**
     * Returns the path to an executable or null if it doesn't exist.
     * @param executable Relative path
     * @return The absolute file path.
     */
    public static File which(String executable) {
        for (String path: PATHS) {
            File file = new File(path, executable);
            if (file.exists())
                return file.getAbsoluteFile();
        }
        return null;
    }

    /**
     * Return the current classpath as a list of absolute paths
     * @return
     */
    public static List<String> getAbsoluteClassPaths() {
        final String[] relativeClassPaths = System.getProperty("java.class.path").split(File.pathSeparator);
        final List<String> absoluteClassPaths = new ArrayList<>(relativeClassPaths.length);
        for (String classPath : relativeClassPaths) {
            File cp = new File(classPath);
            if (cp.exists())
                absoluteClassPaths.add(cp.getAbsolutePath());
        }

        return absoluteClassPaths;
    }
}
