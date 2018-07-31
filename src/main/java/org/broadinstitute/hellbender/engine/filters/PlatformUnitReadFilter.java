package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTag;
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
 * Filter out reads where the the platform unit attribute (PU tag) contains the given string.
 *
 * <p>Note: Matching is done by exact case-sensitive text matching.</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads with matching platform unit attribute")
public final class PlatformUnitReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.BLACK_LISTED_LANES_LONG_NAME,
            doc="Platform unit (PU) to filter out",
            optional=false)
    public Set<String> blackListedLanes = new LinkedHashSet<>();

    // Command line parser requires a no-arg constructor
    public PlatformUnitReadFilter( ) {}

    public PlatformUnitReadFilter( final SAMFileHeader header ) {
        super.setHeader(header);
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

        return ReadUtils.getPlatformUnit(read, samHeader);
    }
}
