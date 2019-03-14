package com.intel.gkl;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.net.URL;
import java.util.HashSet;
import java.util.Set;

/**
 * Loads native libraries from the classpath, usually from a jar file.
 */
public final class NativeLibraryLoader {
    private static final Logger logger = LogManager.getLogger(NativeLibraryLoader.class);
    private static final String USE_LIBRARY_PATH = "USE_LIBRARY_PATH";
    private static final Set<String> loadedLibraries = new HashSet<String>();

    /**
     * Tries to load the native library from the classpath, usually from a jar file. <p>
     *
     * If the USE_LIBRARY_PATH environment variable is defined, the native library will be loaded from the
     * java.library.path instead of the classpath.
     *
     * @param tempDir  directory where the native library is extracted or null to use the system temp directory
     * @param libraryName  name of the shared library without system dependent modifications
     * @return true if the library was loaded successfully, false otherwise
     */
    public static synchronized boolean load(File tempDir, String libraryName) throws Exception {
        if (loadedLibraries.contains(libraryName)) {
            return true;
        }

        final String systemLibraryName = System.mapLibraryName(libraryName);

        // load from the java.library.path
        if (System.getenv(USE_LIBRARY_PATH) != null) {
            final String javaLibraryPath = System.getProperty("java.library.path");
            try {
                logger.error(String.format("OVERRIDE DEFAULT: Loading %s from %s", systemLibraryName, javaLibraryPath));
                logger.error(String.format("LD_LIBRARY_PATH = %s", System.getenv("LD_LIBRARY_PATH")));
                System.loadLibrary(libraryName);
                return true;
            } catch (Exception|Error e) {
                logger.error(String.format("OVERRIDE DEFAULT: Unable to load %s from %s", systemLibraryName, javaLibraryPath));
                throw new RuntimeException(String.format("OVERRIDE DEFAULT: Unable to load %s from %s", systemLibraryName, javaLibraryPath)+e.getCause());
//                return false;
            }
        }

        // load from the java classpath
        final String resourcePath = "native/" +  systemLibraryName;
        final URL inputUrl = NativeLibraryLoader.class.getResource(resourcePath);
        if (inputUrl == null) {
            logger.error("Unable to find native library: " + resourcePath);
            throw new RuntimeException("Unable to find native library: " + resourcePath);
//            return false;
        }
        logger.error(String.format("Loading %s from %s", systemLibraryName, inputUrl.toString()));

        try {
            final File temp = File.createTempFile(FilenameUtils.getBaseName(resourcePath),
                    "." + FilenameUtils.getExtension(resourcePath), tempDir);
            FileUtils.copyURLToFile(inputUrl, temp);
            temp.deleteOnExit();
            logger.error(String.format("Extracting %s to %s", systemLibraryName, temp.getAbsolutePath()));
            System.load(temp.getAbsolutePath());
        } catch (Exception|Error e) {
            logger.error(String.format("Unable to load %s from %s (%s)", systemLibraryName, resourcePath, e.getMessage()));
            throw new RuntimeException(String.format("Unable to load %s from %s (%s)", systemLibraryName, resourcePath, e.getMessage()));
            //            return false;
        }

        loadedLibraries.add(libraryName);
        return true;
    }

}
