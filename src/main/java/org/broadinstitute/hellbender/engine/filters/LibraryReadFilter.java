package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;

/**
 * Keep only reads from the specified library.
 */
public final class LibraryReadFilter implements ReadFilter, CommandLineFilter, Serializable {
    private static final long serialVersionUID = 1L;

    private final static String libraryArgName = "library";

    @Argument(fullName = libraryArgName, shortName = "library", doc="The name of the library to keep", optional=true)
    public String libraryToKeep = null;

    private final SAMFileHeader header;

    public LibraryReadFilter( final SAMFileHeader header ) {
        this.header = header;
    }

    @Override
    public boolean test( final GATKRead read ) {
        final String library = ReadUtils.getLibrary(read, header);
        return library != null && library.equals(libraryToKeep);
    }

    @Override
    public String validate() {
        String message = null;
        if (libraryToKeep == null || libraryToKeep.length() < 1) {
            message = "A value for \"" + libraryArgName + "\" must be provided";
        }
        return message;
    }

}
