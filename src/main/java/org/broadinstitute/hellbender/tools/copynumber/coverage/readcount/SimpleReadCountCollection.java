package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;


import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ioutils.HDF5Utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;

/**
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class SimpleReadCountCollection extends ReadCountCollection<SimpleReadCountData> {

    @Override
    Supplier<RealMatrix> initializeRealMatrix(final File file) {
        if (HDF5Utils.isHDF5File(file)) {
            //read in HDF5 file
        } else {

        }
    }

    @Override
    Supplier<Map<SimpleInterval, SimpleReadCountData>> initializeReadCountDataMap(final File file) {
        return null;
    }

    @Override
    String readSampleName(final File file) {
        return null;
    }
}
