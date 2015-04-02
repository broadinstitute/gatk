package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTag;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.HashSet;
import java.util.Set;

/**
 * Keep reads that do not have blacklisted platform unit tags.
 * Matching is done by exact case-sensitive text matching.
 */
public final class PlatformUnitReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "blackListedLanes", shortName = "blackListedLanes", doc="Keep reads with platform units not on the list", optional=true)
    public Set<String> blackListedLanes = new HashSet<>();

    private final SAMFileHeader header;

    public PlatformUnitReadFilter( final SAMFileHeader header ) {
        this.header = header;
    }

    @Override
    public boolean test( final GATKRead read ) {
        if (blackListedLanes.isEmpty()) {
            return true;
        }
        final String pu_attr = getPlatformUnit(read);
        return pu_attr == null || !blackListedLanes.contains(pu_attr);
    }

    private String getPlatformUnit( final GATKRead read ) {
        final String pu_attr = read.getAttributeAsString(SAMTag.PU.name());
        if ( pu_attr != null ) {
            return pu_attr;
        }

        return ReadUtils.getPlatformUnit(read, header);
    }
}
