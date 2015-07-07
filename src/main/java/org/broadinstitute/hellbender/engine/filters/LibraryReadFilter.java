package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Keep only reads from the specified library.
 */
public final class LibraryReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "library", shortName = "library", doc="The name of the library to keep", optional=false)
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
}
