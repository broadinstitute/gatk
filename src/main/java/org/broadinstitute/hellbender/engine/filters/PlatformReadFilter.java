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
 * Keep only reads where the the Read Group platform attribute (RG:PL tag) contains the given string.
 *
 * <p>Note: Matching is done by case-insensitive substring matching.</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads with matching Read Group platform")
public final class PlatformReadFilter extends ReadFilter implements Serializable{
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.PL_FILTER_NAME_LONG_NAME,
            doc="Platform attribute (PL) to match", optional=false)
    public Set<String> PLFilterNames = new LinkedHashSet<>();

    // Command line parser requires a no-arg constructor
    public PlatformReadFilter() {}

    public PlatformReadFilter( final SAMFileHeader header ) {
        super.setHeader(header);
    }

    @Override
    public boolean test( final GATKRead read ) {
        String platform = ReadUtils.getPlatform(read, samHeader);
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
