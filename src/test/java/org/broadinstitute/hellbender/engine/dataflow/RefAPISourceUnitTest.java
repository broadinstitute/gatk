package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPIMetadata;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Map;

public class RefAPISourceUnitTest extends BaseTest {

    private ReferenceBases queryReferenceAPI( final String referenceName, final SimpleInterval interval ) {
        GenomicsOptions options = PipelineOptionsFactory.create().as(GenomicsOptions.class);
        options.setApiKey(getDataflowTestApiKey());
        options.setProject(getDataflowTestProject());

        final Pipeline p = TestPipeline.create(options); // We don't use GATKTestPipeline because we need specific options.

        Map<String, String> referenceNameToIdTable = RefAPISource.buildReferenceNameToIdTable(p.getOptions(), referenceName);
        RefAPISource refAPISource = RefAPISource.getRefAPISource();
        RefAPIMetadata refAPIMetadata = new RefAPIMetadata(referenceName, referenceNameToIdTable);

        return refAPISource.getReferenceBases(p.getOptions(), refAPIMetadata, interval);
    }

    @Test(groups = "cloud_todo")
    public void testDummy() {
        String referenceName = "EOSt9JOVhp3jkwE";
        final String expected = "AAACAGGTTA";
        // -1 because we're using closed intervals
        SimpleInterval interval = new SimpleInterval("1", 50001, 50001 + expected.length() - 1);
        Logger logger = LogManager.getLogger(RefAPISourceUnitTest.class);

        GenomicsOptions options = PipelineOptionsFactory.create().as(GenomicsOptions.class);
        options.setApiKey(getDataflowTestApiKey());
        options.setProject(getDataflowTestProject());

        final Pipeline p = TestPipeline.create(options); // We don't use GATKTestPipeline because we need specific options.
        logger.info("a");
        Map<String, String> referenceNameToIdTable = RefAPISource.buildReferenceNameToIdTable(p.getOptions(), referenceName);
        logger.info("b");
        RefAPISource refAPISource = RefAPISource.getRefAPISource();
        logger.info("c");
        RefAPIMetadata refAPIMetadata = new RefAPIMetadata(referenceName, referenceNameToIdTable);
        logger.info("d");
        ReferenceBases bases = refAPISource.getReferenceBases(p.getOptions(), refAPIMetadata, interval);
        logger.info("e");

        final String actual = new String(bases.getBases());
        Assert.assertEquals(actual, expected, "Wrong bases returned");
        p.run();

    }

    @Test(groups = "cloud_todo")
    public void testReferenceSourceQuery() {
        final ReferenceBases bases = queryReferenceAPI("EOSt9JOVhp3jkwE", new SimpleInterval("1", 50000, 50009));

        Assert.assertNotNull(bases);
        Assert.assertNotNull(bases.getBases());
        Assert.assertEquals(bases.getBases().length, 10, "Wrong number of bases returned");
        Assert.assertEquals(new String(bases.getBases()), "TAAACAGGTT", "Wrong bases returned");
    }

    @Test(groups = "cloud_todo", expectedExceptions = UserException.class)
    public void testReferenceSourceQueryWithInvalidContig() {
        final ReferenceBases bases = queryReferenceAPI("EOSt9JOVhp3jkwE", new SimpleInterval("FOOCONTIG", 1, 2));
    }

    @Test(groups = "cloud_todo", expectedExceptions = UserException.class)
    public void testReferenceSourceQueryWithInvalidPosition() {
        final ReferenceBases bases = queryReferenceAPI("EOSt9JOVhp3jkwE", new SimpleInterval("1", 1000000000, 2000000000));
    }

    @Test(groups = "cloud_todo", expectedExceptions = IllegalArgumentException.class)
    public void testReferenceSourceQueryWithNullInterval() {
        final ReferenceBases bases = queryReferenceAPI("EOSt9JOVhp3jkwE", null);
    }

}