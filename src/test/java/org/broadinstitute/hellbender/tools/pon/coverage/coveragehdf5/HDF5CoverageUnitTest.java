package org.broadinstitute.hellbender.tools.pon.coverage.coveragehdf5;

import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

public class HDF5CoverageUnitTest extends BaseTest {

    @Test
    public void basicTest() throws IOException {
        final File outputFile = createTempFile("HDF5Coverage", ".cov");
        final List<Target> newTargets = new ArrayList<>();
        newTargets.add(new Target("target1", new SimpleInterval("1", 1000, 2000)));
        newTargets.add(new Target("target2", new SimpleInterval("1", 5000, 6000)));

        final List<String> sampleNames = new ArrayList<>();
        sampleNames.add("SAMPLE1");
        sampleNames.add("SAMPLE2");

        final double[][] newTargetValues = {{2,4}, {10,12}};

        // The output file already exists at this point, since it is a temp file.
        HDF5Coverage.write(outputFile, HDF5File.OpenMode.READ_WRITE, newTargets, newTargetValues, sampleNames);

        final ReadCountCollection rcc = ReadCountCollectionUtils.parseHdf5AsDouble(outputFile);
        Assert.assertEquals(rcc.columnNames(), sampleNames);
        Assert.assertEquals(rcc.targets(), newTargets);
        Assert.assertFalse(rcc.targets() == newTargets);
        Assert.assertFalse(rcc.columnNames() == sampleNames);

        for (int j = 0; j < newTargetValues.length; j ++) {
            final int j0 = j;
            Assert.assertTrue(IntStream.range(0, newTargetValues[j0].length).allMatch(i -> newTargetValues[j0][i] == rcc.counts().getData()[j0][i]));
        }

        Assert.assertEquals(rcc.counts().getData().length, newTargetValues.length);
        Assert.assertEquals(rcc.counts().getData()[0].length, newTargetValues[0].length);
        Assert.assertFalse(rcc.counts().getData() == newTargetValues);
    }
}
