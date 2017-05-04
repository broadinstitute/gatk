package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Filters out reads marked as duplicates.
 */
public final class MarkDuplicatesReadFilter extends ReadFilter implements Serializable {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "Only filter optical duplicates",
            fullName = "opticalOnly",
            optional = true)
    public boolean filterOpticalOnly = false;

    public MarkDuplicatesReadFilter() { }

    public MarkDuplicatesReadFilter(final boolean filterOpticalOnly) {
        this.filterOpticalOnly = filterOpticalOnly;
    }

    @Override
    public boolean test(final GATKRead read) {
        return !((!filterOpticalOnly && read.isDuplicate()) || (read.hasAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)
                && read.getAttributeAsInteger(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME) != 0));
    }
}
