package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Filters out reads marked as duplicates.
 */
public final class NotOpticalDuplicateReadFilter extends ReadFilter implements Serializable {

    private static final long serialVersionUID = 1L;

    public NotOpticalDuplicateReadFilter() { }

    @Override
    public boolean test(final GATKRead read) {
        final Integer odTag = read.getAttributeAsInteger(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME);
        return odTag == null || odTag.intValue() == 0;
    }
}
