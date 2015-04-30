package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;


public class ReferenceSourceUnitTest {

    final String API_KEY = System.getenv("GOOGLE_API_KEY");
    final String TEST_PROJECT = System.getenv("TEST_PROJECT");

    @Test
    public void testReferenceSourceQuery() {
        final String referenceName = "EOSt9JOVhp3jkwE"; //"GRCh37";
        final Pipeline p = TestPipeline.create();

        GenomicsOptions options = PipelineOptionsFactory.create().as(GenomicsOptions.class);
        options.setApiKey(API_KEY);
        options.setProject(TEST_PROJECT);

        final ReferenceSource refSource = new ReferenceSource(referenceName, options);
        final ReferenceBases bases = refSource.getReferenceBases(new SimpleInterval("1", 50000, 50010));

        Assert.assertNotNull(bases);
        Assert.assertNotNull(bases.getBases());
        Assert.assertEquals(bases.getBases().length, 10, "Wrong number of bases returned");
        Assert.assertEquals(new String(bases.getBases()), "TAAACAGGTT", "Wrong bases returned");
    }

}
