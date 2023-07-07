package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Set;

/**
 * Keep only reads with one of these readnames.
 *
 * Should you need to provide many readnames to keep, you can provide a newline separated .list file of readnames to keep.
 *
 * <p>Matching is done by case-sensitive exact match.</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads with this read name")
public final class ReadNameReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, doc="Keep only reads with this read name", optional=false)
    public Set<String> readNames = new LinkedHashSet<>(0);

    @Override
    public boolean test( final GATKRead read ) {
        return read.getName() != null && readNames.contains(read.getName());
    }

}
