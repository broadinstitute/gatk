package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class RefAPISourceUnitTest extends BaseTest {

    private ReferenceBases queryReferenceAPI( final String referenceName, final SimpleInterval interval ) {
        GenomicsOptions options = PipelineOptionsFactory.create().as(GenomicsOptions.class);
        options.setApiKey(getDataflowTestApiKey());
        options.setProject(getDataflowTestProject());

        final Pipeline p = TestPipeline.create(options); // We don't use GATKTestPipeline because we need specific options.

        RefAPISource refAPISource = new RefAPISource(referenceName, p.getOptions());
        return refAPISource.getReferenceBases(interval, p.getOptions());
    }

    @Test(groups = "cloud_todo")
    public void testDummy() {
        String referenceName = "EOSt9JOVhp3jkwE";
        SimpleInterval interval = new SimpleInterval("1", 50001, 10050000);
        Logger logger = LogManager.getLogger(RefAPISourceUnitTest.class);

        GenomicsOptions options = PipelineOptionsFactory.create().as(GenomicsOptions.class);
        options.setApiKey(getDataflowTestApiKey());
        options.setProject(getDataflowTestProject());

        final Pipeline p = TestPipeline.create(options); // We don't use GATKTestPipeline because we need specific options.
        logger.info("a");
        RefAPISource refAPISource = new RefAPISource(referenceName, p.getOptions());
        logger.info("b");
        ReferenceBases bases = refAPISource.getReferenceBases(interval, p.getOptions());
        logger.info("c");

        Assert.assertEquals(new String(bases.getBases()), "AAACAGGTTA", "Wrong bases returned");
        p.run();

    }

    @Test(groups = "cloud")
    public void testReferenceSourceQuery() {
        final ReferenceBases bases = queryReferenceAPI("EOSt9JOVhp3jkwE", new SimpleInterval("1", 50000, 50009));

        Assert.assertNotNull(bases);
        Assert.assertNotNull(bases.getBases());
        Assert.assertEquals(bases.getBases().length, 10, "Wrong number of bases returned");
        Assert.assertEquals(new String(bases.getBases()), "TAAACAGGTT", "Wrong bases returned");
    }

    @Test(groups = "cloud", expectedExceptions = UserException.class)
    public void testReferenceSourceQueryWithInvalidContig() {
        final ReferenceBases bases = queryReferenceAPI("EOSt9JOVhp3jkwE", new SimpleInterval("FOOCONTIG", 1, 2));
    }

    @Test(groups = "cloud", expectedExceptions = UserException.class)
    public void testReferenceSourceQueryWithInvalidPosition() {
        final ReferenceBases bases = queryReferenceAPI("EOSt9JOVhp3jkwE", new SimpleInterval("1", 1000000000, 2000000000));
    }

    @Test(groups = "cloud", expectedExceptions = IllegalArgumentException.class)
    public void testReferenceSourceQueryWithNullInterval() {
        final ReferenceBases bases = queryReferenceAPI("EOSt9JOVhp3jkwE", null);
    }

}