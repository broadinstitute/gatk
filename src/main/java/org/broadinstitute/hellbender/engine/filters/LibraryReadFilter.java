package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.LinkedHashSet;
import java.util.Set;

/**
 * Keep only reads from the specified library.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads from the specified library")
public final class LibraryReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.LIBRARY_NAME, shortName = ReadFilterArgumentDefinitions.LIBRARY_NAME, doc="Name of the library to keep", optional=false)
    public Set<String> libraryToKeep = new LinkedHashSet<>();

    // Command line parser requires a no-arg constructor
    public LibraryReadFilter() { }
    public LibraryReadFilter( final SAMFileHeader header ) {
        super.setHeader(header);
    }

    @Override
    public boolean test( final GATKRead read ) {
        final String library = ReadUtils.getLibrary(read, samHeader);
        if (library == null) {
            return false;
        }

        return libraryToKeep.contains(library);
    }
}
