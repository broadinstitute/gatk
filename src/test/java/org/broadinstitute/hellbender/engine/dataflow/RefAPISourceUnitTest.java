package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;
import com.google.api.services.genomics.model.Reference;

public class RefAPISourceUnitTest extends BaseTest {

    private ReferenceBases queryReferenceAPI( final String referenceName, final SimpleInterval interval ) {
        GenomicsOptions options = PipelineOptionsFactory.create().as(GenomicsOptions.class);
        options.setApiKey(getDataflowTestApiKey());
        options.setProject(getDataflowTestProject());

        final Pipeline p = TestPipeline.create(options); // We don't use GATKTestPipeline because we need specific options.

        RefAPISource refAPISource = new RefAPISource(p.getOptions(), RefAPISource.URL_PREFIX + referenceName);
        return refAPISource.getReferenceBases(p.getOptions(), interval);
    }

    @DataProvider(name = "sortData")
    public Object[][] createSortData() {
        return new String[][] {
            // numerical order, not alphabetic.
            { "1,10,2",
                "1,2,10" },
            // x sorted at the right place
            { "y,x,1,10,2",
                "1,2,10,x,y" },
            // mitochondrial at the end
            { "mt,11,12,13,14,x",
                "11,12,13,14,x,mt" },
            // order works also for 'chr'
            { "chr1,chr10,chr2",
                "chr1,chr2,chr10" },
            // order works also for 'ch'
            { "ch1,ch10,ch2,x",
                "ch1,ch2,ch10,x" },
        };
    }

    @Test(dataProvider="sortData")
    public void testSequenceDictionarySorting(String inputs, String outputs) {
        final String[] input = inputs.split(",");
        final String[] expected = outputs.split(",");
        final RefAPISource ref = new RefAPISource(createDummyReferenceMap(input));
        final SAMSequenceDictionary seq = ref.getReferenceSequenceDictionary(null);
        checkSequenceDictionary(seq, expected);
    }



    @Test(groups = "cloud")
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
        RefAPISource refAPISource = new RefAPISource(p.getOptions(), RefAPISource.URL_PREFIX + referenceName);
        ReferenceBases bases = refAPISource.getReferenceBases(p.getOptions(), interval);
        final String actual = new String(bases.getBases());
        Assert.assertEquals(actual, expected, "Wrong bases returned");
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


    private Map<String, Reference> createDummyReferenceMap(String[] contig) {
        HashMap<String,Reference> fake = new HashMap<>();
        for (String s : contig) {
            Reference r = new Reference();
            String id = "id-"+s;
            r.setName(s);
            r.setId(id);
            r.setLength(100L);
            fake.put(s, r);
        }
        return fake;
    }

    private void checkSequenceDictionary(SAMSequenceDictionary seq, String[] contig) {
        int i=0;
        for (String s: contig) {
            Assert.assertEquals(seq.getSequence(i).getSequenceName(), s);
            i++;
        }
    }

}
