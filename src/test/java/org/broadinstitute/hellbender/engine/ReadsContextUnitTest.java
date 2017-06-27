package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.ReadNameReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Iterator;

public final class ReadsContextUnitTest extends BaseTest {

    @DataProvider(name = "EmptyReadsContextDataProvider")
    public Object[][] getEmptyReadsContextData() {
        // Default-constructed ReadsContexts and ReadsContexts constructed from null ReadsDataSources/intervals
        // should behave as empty context objects.
        return new Object[][] {
                { new ReadsContext() },
                { new ReadsContext(null, null) },
                { new ReadsContext(null, new SimpleInterval("1", 1, 1) ) },
                { new ReadsContext(new ReadsDataSource(IOUtils.getPath(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam")), null) }
        };
    }

    @Test(dataProvider = "EmptyReadsContextDataProvider")
    public void testEmptyReadsContext( final ReadsContext readsContext ) {
        Assert.assertFalse(readsContext.hasBackingDataSource() && readsContext.getInterval() != null,
                           "Empty ReadsContext reports having both a backing data source and an interval");
        Assert.assertFalse(readsContext.iterator().hasNext(), "Empty ReadsContext should have returned an empty iterator from iterator()");
    }

    @DataProvider(name = "ReadsContextFilteringDataProvider")
    public Object[][] getValidReadsContextData() {
        // Default-constructed ReadsContexts and ReadsContexts constructed from null ReadsDataSources/intervals
        // should behave as empty context objects.
        ReadNameReadFilter readNameFilter = new ReadNameReadFilter();
        readNameFilter.readName = "d";
        return new Object[][]{
                {new ReadsContext(
                        new ReadsDataSource(IOUtils.getPath(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam")),
                        new SimpleInterval("1", 200, 210), // query over small interval, with no read filter
                        ReadFilterLibrary.ALLOW_ALL_READS), new String[] { "a", "b", "c" },
                },
                {new ReadsContext(
                        new ReadsDataSource(IOUtils.getPath(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam")),
                        new SimpleInterval("1", 200, 1000), // query over larger interval with readNameFilter on "d"
                        readNameFilter), new String[] { "d" }
                }
        };
    }

    @Test(dataProvider = "ReadsContextFilteringDataProvider")
    public void testReadsContextFiltering(final ReadsContext readsContext, String[] expectedReads ) {
        Assert.assertTrue(readsContext.hasBackingDataSource() && readsContext.getInterval() != null,
                "Valid ReadsContext reports having no backing data source or interval");
        Iterator<GATKRead> it = readsContext.iterator();
        while (it.hasNext()) {
            Assert.assertTrue(Arrays.asList(expectedReads).contains(it.next().getName()));
        }
    }

}