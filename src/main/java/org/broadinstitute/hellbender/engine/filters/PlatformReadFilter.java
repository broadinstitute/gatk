package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Keep only reads that match th PL attribute.
 * Matching is done by case-insensitive substring matching
 * (checking if the read's platform tag contains the given string).
 */
public final class PlatformReadFilter implements ReadFilter, CommandLineFilter, Serializable {
    private static final long serialVersionUID = 1L;

    private final static String platformArgName = "PLFilterName";

    @Argument(fullName = platformArgName, shortName = "PLFilterName", doc="Keep reads with RG:PL attribute containing this string", optional=true)
    public Set<String> PLFilterNames = new LinkedHashSet<>();

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

    @Override
    public String validate() {
        String message = null;
        if (PLFilterNames.size() <= 0) {
            message = "requires one or more values for \"" + platformArgName + "\"";
        }

        return message;
    }

}
