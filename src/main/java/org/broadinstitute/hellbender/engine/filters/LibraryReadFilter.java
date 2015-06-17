package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;

/**
 * Keep only reads from the specified library.
 */
public final class LibraryReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "library", shortName = "library", doc="The name of the library to keep", optional=false)
    public String libraryToKeep = null;

    @Override
    public boolean test(final SAMRecord read) {
        final SAMReadGroupRecord readGroup = read.getReadGroup();
        return readGroup != null && readGroup.getLibrary() != null && readGroup.getLibrary().equals(libraryToKeep);
    }
}
