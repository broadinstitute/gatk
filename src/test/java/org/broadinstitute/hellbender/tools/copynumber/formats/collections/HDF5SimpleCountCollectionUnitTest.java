package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class HDF5SimpleCountCollectionUnitTest extends GATKBaseTest {
    @Test
    public void basicTest() {
        final File outputFile = createTempFile("HDF5ReadCountCollection", ".hdf5");
        final SampleLocatableMetadata metadata = new SimpleSampleLocatableMetadata(
                "test-sample",
                new SAMSequenceDictionary(Collections.singletonList(
                        new SAMSequenceRecord("1", 10000))));
        final List<SimpleInterval> intervals = new ArrayList<>();
        intervals.add(new SimpleInterval("1", 1000, 2000));
        intervals.add(new SimpleInterval("1", 5000, 6000));
        final double[] counts = {2, 10};
        HDF5SimpleCountCollection.write(outputFile, metadata, intervals, counts);

        final HDF5SimpleCountCollection rcc = new HDF5SimpleCountCollection(new HDF5File(outputFile));
        Assert.assertEquals(rcc.getMetadata(), metadata);
        Assert.assertEquals(rcc.getIntervals(), intervals);
        Assert.assertEquals(rcc.getCounts().getRowDimension(), 1);
        Assert.assertFalse(rcc.getCounts().getRow(0) == counts);
    }
}
