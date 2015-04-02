package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Set;

/**
 * Keep only reads that match th PL attribute.
 * Matching is done by case-insensitive substring matching
 * (checking if the read's platform tag contains the given string).
 */
public final class PlatformReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "PLFilterName", shortName = "PLFilterName", doc="Keep reads with RG:PL attribute containing this string", optional=true)
    public Set<String> PLFilterNames;

    private final SAMFileHeader header;

    public PlatformReadFilter( final SAMFileHeader header ) {
        this.header = header;
    }

    @Override
    public boolean test( final GATKRead read ) {
        String platform = ReadUtils.getPlatform(read, header);
        if ( platform == null ) {
            return false;
        }
        platform = platform.toUpperCase();

        for ( final String name : PLFilterNames ) {
            if ( platform.contains(name.toUpperCase()) ) {
                return true;
            }
        }
        return false;
    }
}
