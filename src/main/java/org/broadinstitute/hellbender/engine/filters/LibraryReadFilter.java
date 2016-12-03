package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;

/**
 * Keep only reads from the specified library.
 */
public final class LibraryReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "library", shortName = "library", doc="The name of the library to keep", optional=false)
    public String libraryToKeep = null;

    // Command line parser requires a no-arg constructor
    public LibraryReadFilter() { }
    public LibraryReadFilter( final SAMFileHeader header ) {
        super.setHeader(header);
    }

    @Override
    public boolean test( final GATKRead read ) {
        final String library = ReadUtils.getLibrary(read, samHeader);
        return library != null && library.equals(libraryToKeep);
    }
}
