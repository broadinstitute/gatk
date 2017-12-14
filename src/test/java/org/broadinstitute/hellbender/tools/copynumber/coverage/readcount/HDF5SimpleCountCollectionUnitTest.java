package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public final class HDF5SimpleCountCollectionUnitTest extends BaseTest {
    @Test
    public void basicTest() {
        final File outputFile = createTempFile("HDF5ReadCountCollection", ".hdf5");
        final String sampleName = "SAMPLE1";
        final List<SimpleInterval> intervals = new ArrayList<>();
        intervals.add(new SimpleInterval("1", 1000, 2000));
        intervals.add(new SimpleInterval("1", 5000, 6000));
        final double[] counts = {2, 10};

        // The output file already exists at this point, since it is a temp file.
        HDF5SimpleCountCollection.write(outputFile, sampleName, intervals, counts);

        final HDF5SimpleCountCollection rcc = new HDF5SimpleCountCollection(new HDF5File(outputFile));

        Assert.assertEquals(rcc.getSampleName(), sampleName);

        Assert.assertEquals(rcc.getIntervals(), intervals);
        Assert.assertFalse(rcc.getIntervals() == intervals);

        Assert.assertEquals(rcc.getCounts().getRowDimension(), 1);
        Assert.assertFalse(rcc.getCounts().getRow(0) == counts);
    }
}
