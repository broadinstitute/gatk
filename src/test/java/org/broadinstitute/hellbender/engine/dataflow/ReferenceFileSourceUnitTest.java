package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Map;


public class ReferenceFileSourceUnitTest extends BaseTest {

    @Test
    public void testReferenceSourceQuery() throws IOException {
        final String referencePath = hg19MiniReference;

        final ReferenceFileSource refSource = new ReferenceFileSource(referencePath);
        final ReferenceBases bases = refSource.getReferenceBases(PipelineOptionsFactory.create(),
                new SimpleInterval("2", 10001, 10010));

        Assert.assertNotNull(bases);
        Assert.assertEquals(bases.getBases().length, 10, "Wrong number of bases returned");
        Assert.assertEquals(new String(bases.getBases()), "CGTATCCCAC", "Wrong bases returned");
    }

    @Test
    public void testGetReferenceBasesByContig() throws IOException {
        final String referencePath = hg19MiniReference;

        final ReferenceFileSource refSource = new ReferenceFileSource(referencePath);
        Map<String, ReferenceBases> ref = refSource.getAllReferenceBases();

        Assert.assertNotNull(ref);
        Assert.assertEquals(ref.size(), 4);
        ReferenceBases bases = ref.get("2");
        Assert.assertNotNull(bases);
        Assert.assertEquals(bases.getBases().length, 16000, "Wrong number of bases returned");
    }

}
