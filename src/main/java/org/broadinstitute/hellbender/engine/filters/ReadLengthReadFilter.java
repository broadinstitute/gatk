package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep only reads whose length is &ge; min value and &le; max value.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads whose length is within a certain range")
public final class ReadLengthReadFilter extends ReadFilter implements Serializable{
    private static final long serialVersionUID = 1L;

    private static final String maxLengthArgName = "maxReadLength";

    @Argument(fullName = maxLengthArgName,
            shortName = maxLengthArgName,
            doc="Keep only reads with length at most equal to the specified value",
            optional=false)
    public Integer maxReadLength;

    private static final String minLengthArg = "minReadLength";
    @Argument(fullName = minLengthArg,
            shortName = minLengthArg,
            doc="Keep only reads with length at least equal to the specified value",
            optional=true)
    public int minReadLength = 1;

    // Command line parser requires a no-arg constructor
    public ReadLengthReadFilter() {}

    public ReadLengthReadFilter( final int minLength, final int maxLength ) {
        this.minReadLength = minLength;
        this.maxReadLength = maxLength;
    }

    @Override
    public boolean test( final GATKRead read ) {
        return read.getLength() >= minReadLength && read.getLength() <= maxReadLength;
    }

}
