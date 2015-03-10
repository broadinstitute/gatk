package org.broadinstitute.hellbender.utils.R;

import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;

/**
 * Libraries embedded in the StingUtils package.
 */
public enum RScriptLibrary {
    GSALIB("gsalib");

    private final String name;

    private RScriptLibrary(String name) {
        this.name = name;
    }

    public String getLibraryName() {
        return this.name;
    }

    public String getResourcePath() {
        return name + ".tar.gz";
    }

    /**
     * Writes the library source code to a temporary tar.gz file and returns the path.
     * @return The path to the library source code. The caller must delete the code when done.
     */
    public File writeTemp() {
        return IOUtils.writeTempResource(new Resource(getResourcePath(), RScriptLibrary.class));
    }

    public File writeLibrary(File tempDir) {
        File libraryFile = new File(tempDir, getLibraryName());
        IOUtils.writeResource(new Resource(getResourcePath(), RScriptLibrary.class), libraryFile);
        return libraryFile;
    }
}
