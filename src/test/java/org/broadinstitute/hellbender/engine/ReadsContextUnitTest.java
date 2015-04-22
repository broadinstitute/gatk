package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class ReadsContextUnitTest extends BaseTest {

    @DataProvider(name = "EmptyReadsContextDataProvider")
    public Object[][] getEmptyReadsContextData() {
        // Default-constructed ReadsContexts and ReadsContexts constructed from null ReadsDataSources/intervals
        // should behave as empty context objects.
        return new Object[][] {
                { new ReadsContext() },
                { new ReadsContext(null, null) },
                { new ReadsContext(null, new SimpleInterval("1", 1, 1) ) },
                { new ReadsContext(new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam")), null) }
        };
    }

    @Test(dataProvider = "EmptyReadsContextDataProvider")
    public void testEmptyReadsContext( final ReadsContext readsContext ) {
        Assert.assertFalse(readsContext.hasBackingDataSource() && readsContext.getInterval() != null,
                           "Empty ReadsContext reports having both a backing data source and an interval");
        Assert.assertFalse(readsContext.iterator().hasNext(), "Empty ReadsContext should have returned an empty iterator from iterator()");
    }

}
