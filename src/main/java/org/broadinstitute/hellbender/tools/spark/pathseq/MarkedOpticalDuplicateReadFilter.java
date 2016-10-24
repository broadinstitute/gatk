package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

public final class MarkedOpticalDuplicateReadFilter extends ReadFilter implements Serializable {

    private static final long serialVersionUID = 1L;

    //Filters out optical duplicate reads (marked with the OD tag)
    @Override
    public boolean test( final GATKRead read ) {
        if (read.hasAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)
                && read.getAttributeAsInteger(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME) != 0) {
            return false;
        }
        return true;
    }
}
