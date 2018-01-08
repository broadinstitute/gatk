package org.broadinstitute.hellbender.engine.datasources;

import com.google.api.services.genomics.model.Reference;
import com.google.cloud.genomics.utils.GenomicsFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.LinkedHashMap;
import java.util.Map;

import static org.broadinstitute.hellbender.engine.datasources.ReferenceAPISource.*;

public class ReferenceAPISourceUnitTest extends GATKBaseTest {

    private ReferenceBases queryReferenceAPI(final String referenceName, final SimpleInterval interval, int pageSize ) {
        ReferenceAPISource refAPISource = makeReferenceAPISource(referenceName);
        return refAPISource.getReferenceBases(interval, pageSize);
    }

    private ReferenceBases queryReferenceAPI( final String referenceName, final SimpleInterval interval ) {
        ReferenceAPISource refAPISource = makeReferenceAPISource(referenceName);
        return refAPISource.getReferenceBases(interval);
    }

    private ReferenceAPISource makeReferenceAPISource(String referenceName) {
        return new ReferenceAPISource(ReferenceAPISource.URL_PREFIX + referenceName);
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
        final ReferenceAPISource ref = new ReferenceAPISource(createDummyReferenceMap(input));
        final SAMSequenceDictionary seq = ref.getReferenceSequenceDictionary(null);
        checkSequenceDictionary(seq, expected);
    }

    @DataProvider(name = "assemblyIDData")
    public Object[][] assemblyIDData() {
        return new String[][] {
                { HS37D5_ASSEMBLY_ID, HS37D5_REF_ID },
                { GRCH37_LITE_ASSEMBLY_ID, GRCH37_LITE_REF_ID},
                { GRCH38_ASSEMBLY_ID, GRCH38_REF_ID},
                { HG19_ASSEMBLY_ID, HG19_REF_ID},
        };
    }

    @Test(groups = "bucket", dataProvider = "assemblyIDData")
    public void testCreateByAssemblyID(final String assemblyID, final String refID) {
        final ReferenceAPISource apiSourceByAssemblyName = ReferenceAPISource.fromReferenceSetAssemblyID(assemblyID);
        final ReferenceAPISource apiSourceByID = makeReferenceAPISource(refID);
        Assert.assertEquals(apiSourceByAssemblyName.getReferenceMap(), apiSourceByID.getReferenceMap());
    }

    @DataProvider(name = "assemblyIDDataMultiple")
    public Object[][] assemblyIDDataMultiple() {
        return new String[][] {
                { GRCH37_ASSEMBLY_ID, GRCH37_REF_ID},
        };
    }

    @Test(groups = "bucket", dataProvider = "assemblyIDDataMultiple", expectedExceptions = UserException.MultipleReferenceSets.class)
    public void testCreateByAssemblyIDMultipleReferenceSets(final String assemblyID, final String refID) throws Exception {
        final ReferenceAPISource apiSourceByAssemblyName = ReferenceAPISource.fromReferenceSetAssemblyID( assemblyID);
        final ReferenceAPISource apiSourceByID = makeReferenceAPISource(refID);
        Assert.assertEquals(apiSourceByAssemblyName.getReferenceMap(), apiSourceByID.getReferenceMap());
    }

    @Test(groups = "bucket")
    public void testDummy() {
        String referenceName = HS37D5_REF_ID;
        final String expected = "AAACAGGTTA";
        // -1 because we're using closed intervals
        SimpleInterval interval = new SimpleInterval("1", 50001, 50001 + expected.length() - 1);
        Logger logger = LogManager.getLogger(ReferenceAPISourceUnitTest.class);

        ReferenceAPISource refAPISource = makeReferenceAPISource(referenceName);
        ReferenceBases bases = refAPISource.getReferenceBases( interval);
        final String actual = new String(bases.getBases());
        Assert.assertEquals(actual, expected, "Wrong bases returned");
    }

    @Test(groups = "bucket")
    public void testReferenceSourceQuery() {
        final ReferenceBases bases = queryReferenceAPI(HS37D5_REF_ID, new SimpleInterval("1", 50000, 50009));

        Assert.assertNotNull(bases);
        Assert.assertNotNull(bases.getBases());
        Assert.assertEquals(bases.getBases().length, 10, "Wrong number of bases returned");
        Assert.assertEquals(new String(bases.getBases()), "TAAACAGGTT", "Wrong bases returned");
    }

    @Test(groups = "bucket")
    public void testReferenceSourceMultiPageQuery() {
        final int mio = 1_000_000;
        final ReferenceBases bases1 = queryReferenceAPI(HS37D5_REF_ID, new SimpleInterval("1", 50000, 50000 + mio + 50));
        final ReferenceBases bases2 = queryReferenceAPI(HS37D5_REF_ID, new SimpleInterval("1", 50025, 50025 + mio + 50));

        Assert.assertNotNull(bases1);
        Assert.assertNotNull(bases1.getBases());
        Assert.assertNotNull(bases2);
        Assert.assertNotNull(bases2.getBases());
        // those SimpleIntervals include the end, hence +1
        Assert.assertEquals(bases1.getBases().length, mio + 50 + 1, "Wrong number of bases returned");
        Assert.assertEquals(bases2.getBases().length, mio + 50 + 1, "Wrong number of bases returned");

        // grab some bases around the seam
        ReferenceBases seam1 = bases1.getSubset(new SimpleInterval("1", 50000 + mio - 100, 50000 + mio + 50));
        ReferenceBases seam2 = bases2.getSubset(new SimpleInterval("1", 50000 + mio - 100, 50000 + mio + 50));

        Assert.assertEquals(seam1.getBases(), seam2.getBases(), "seam doesn't match (paging bug?)");

    }

    @Test(groups = "bucket")
    public void testReferenceSourceMultiSmallPagesQuery() {
        int pageSize = 300;
        // not a multiple of pageSize (testing the fetching of a partial page)
        final ReferenceBases bases1 = queryReferenceAPI(HS37D5_REF_ID, new SimpleInterval("1", 50000, 51000), pageSize);
        // multiple of pageSize (testing ending on an exact page boundary)
        final ReferenceBases bases2 = queryReferenceAPI(HS37D5_REF_ID, new SimpleInterval("1", 50025, 50924), pageSize);

        Assert.assertNotNull(bases1);
        Assert.assertNotNull(bases1.getBases());
        Assert.assertNotNull(bases2);
        Assert.assertNotNull(bases2.getBases());
        // those SimpleIntervals include the end, hence +1
        Assert.assertEquals(bases1.getBases().length, 1001, "Wrong number of bases returned");
        Assert.assertEquals(bases2.getBases().length, 900, "Wrong number of bases returned");

        // grab some bases they should have in common
        ReferenceBases seam1 = bases1.getSubset(new SimpleInterval("1", 50025, 50902));
        ReferenceBases seam2 = bases2.getSubset(new SimpleInterval("1", 50025, 50902));

        Assert.assertEquals(seam1.getBases(), seam2.getBases(), "seam doesn't match (paging bug?)");
    }

    @Test(groups = "bucket", expectedExceptions = UserException.ReferenceAPIReturnedUnexpectedNumberOfBytes.class)
    public void testOffContig() throws Exception {
        //Test the case of a query that starts before the end of the contig and runs off
        final int b37Chr1Len = 249250621;
        final int pageSize = 300;

        final int beforeEnd = 11;
        final int pastEnd = 12;
        final int start = b37Chr1Len - beforeEnd;
        final int end   = b37Chr1Len + pastEnd;
        final ReferenceBases bases = queryReferenceAPI(ReferenceAPISource.HS37D5_REF_ID, new SimpleInterval("1", start, end), pageSize);
    }

    @Test(groups = "bucket")
    public void testReferenceSourceVaryingPageSizeQuery() {

        SimpleInterval interval = new SimpleInterval("1", 50000, 50050);
        final ReferenceBases bases1 = queryReferenceAPI(HS37D5_REF_ID, interval);
        final ReferenceBases bases2 = queryReferenceAPI(HS37D5_REF_ID, interval, 10);

        Assert.assertNotNull(bases1);
        Assert.assertNotNull(bases1.getBases());
        Assert.assertNotNull(bases2);
        Assert.assertNotNull(bases2.getBases());
        Assert.assertEquals(bases1.getBases(), bases2.getBases(), "bases should match despite different paging size");
    }

    @Test(groups = "bucket", expectedExceptions = UserException.class)
    public void testReferenceSourceQueryWithInvalidContig() {
        final ReferenceBases bases = queryReferenceAPI(HS37D5_REF_ID, new SimpleInterval("FOOCONTIG", 1, 2));
    }

    @Test(groups = "bucket", expectedExceptions = UserException.class)
    public void testReferenceSourceQueryWithInvalidPosition() {
        final ReferenceBases bases = queryReferenceAPI(HS37D5_REF_ID, new SimpleInterval("1", 1000000000, 2000000000));
    }

    @Test(groups = "bucket", expectedExceptions = IllegalArgumentException.class)
    public void testReferenceSourceQueryWithNullInterval() {
        final ReferenceBases bases = queryReferenceAPI(HS37D5_REF_ID, null);
    }

    private Map<String, Reference> createDummyReferenceMap(String[] contig) {
        Map<String,Reference> fake = new LinkedHashMap<>();
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
