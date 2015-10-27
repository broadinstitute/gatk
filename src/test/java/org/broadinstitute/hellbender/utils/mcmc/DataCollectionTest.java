package org.broadinstitute.hellbender.utils.mcmc;

import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

/**
 * Unit tests for {@link DataCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class DataCollectionTest {
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDuplicateDatasetNamesException() {
        new DataCollection(Arrays.asList(new Data<>("data", Collections.singletonList(1.)),
                new Data<>("data", Collections.singletonList(1.))));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddDuplicateDatasetNamesException() {
        final DataCollection dataCollection =
                new DataCollection(Collections.singletonList(new Data<>("data", Collections.singletonList(1.))));
        dataCollection.add(new Data<>("data", Collections.singletonList(1.)));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetBadDatasetNameException() {
        final DataCollection dataCollection =
                new DataCollection(Collections.singletonList(new Data<>("data1", Collections.singletonList(1.))));
        dataCollection.get("data2", Double.class);
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testGetBadDatasetTypeException() {
        final DataCollection dataCollection =
                new DataCollection(Collections.singletonList(new Data<>("data1", Collections.singletonList(1))));
        dataCollection.get("data1", Double.class);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTypeUnsafeGetBadDatasetNameException() {
        final DataCollection dataCollection =
                new DataCollection(Collections.singletonList(new Data<>("data1", Collections.singletonList(1.))));
        dataCollection.get("data2");
    }
}