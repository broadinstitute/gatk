package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public final class ReferenceDataSourceUnitTest extends GATKBaseTest {

    private static final Path TEST_REFERENCE = IOUtils.getPath(hg19MiniReference);
    private static final Path TEST_REFERENCE_BGZ = IOUtils.getPath(hg19MiniReference + ".gz");

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentReference() {
        new ReferenceFileSource(GATKBaseTest.getSafeNonExistentPath("nonexistent.fasta"));
    }

    @Test(expectedExceptions = UserException.MissingReferenceFaiFile.class)
    public void testReferenceWithMissingFaiFile() {
        ReferenceDataSource refDataSource = new ReferenceFileSource(IOUtils.getPath(publicTestDir + "fastaWithoutFai.fasta"));
    }

    @Test(expectedExceptions = UserException.MissingReferenceDictFile.class)
    public void testReferenceWithMissingDictFile() {
        ReferenceDataSource refDataSource = new ReferenceFileSource(IOUtils.getPath(publicTestDir + "fastaWithoutDict.fasta"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullReference() {
        ReferenceDataSource refDataSource = new ReferenceFileSource((Path)null);
    }

    @Test
    public void testGetSequenceDictionary() {
        try (ReferenceDataSource refDataSource = new ReferenceFileSource(TEST_REFERENCE)) {
            SAMSequenceDictionary sequenceDictionary = refDataSource.getSequenceDictionary();
            Assert.assertEquals(sequenceDictionary.size(), 4, "Wrong number of sequences in sequence dictionary returned from refDataSource.getSequenceDictionary()");
            for ( String contig : Arrays.asList("1", "2", "3", "4") ) {
                Assert.assertNotNull(sequenceDictionary.getSequence(contig), "Sequence dictionary returned from refDataSource.getSequenceDictionary() lacks expected contig " + contig);
            }
        }
    }

    @DataProvider(name = "ReferenceIntervalDataProvider")
    public Object[][] getReferenceIntervals() {
        return new Object[][] {
            { TEST_REFERENCE, new SimpleInterval("1", 1, 3), "NNN" },
            { TEST_REFERENCE, new SimpleInterval("1", 11041, 11045), "GCAAA" },
            { TEST_REFERENCE, new SimpleInterval("1", 11210, 11220), "CGGTGCTGTGC" },
            { TEST_REFERENCE, new SimpleInterval("2", 9995, 10005), "NNNNNNCGTAT" },
            { TEST_REFERENCE, new SimpleInterval("2", 10001, 10080), "CGTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCAC" },
            { TEST_REFERENCE, new SimpleInterval("2", 10005, 10084), "TCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCACACCC" },
            { TEST_REFERENCE, new SimpleInterval("2", 15995, 16000), "TGTCAG" },

            { TEST_REFERENCE_BGZ, new SimpleInterval("1", 1, 3), "NNN" },
            { TEST_REFERENCE_BGZ, new SimpleInterval("1", 11041, 11045), "GCAAA" },
            { TEST_REFERENCE_BGZ, new SimpleInterval("1", 11210, 11220), "CGGTGCTGTGC" },
            { TEST_REFERENCE_BGZ, new SimpleInterval("2", 9995, 10005), "NNNNNNCGTAT" },
            { TEST_REFERENCE_BGZ, new SimpleInterval("2", 10001, 10080), "CGTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCAC" },
            { TEST_REFERENCE_BGZ, new SimpleInterval("2", 10005, 10084), "TCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCACACCC" },
            { TEST_REFERENCE_BGZ, new SimpleInterval("2", 15995, 16000), "TGTCAG" }
        };
    }


    @Test(dataProvider = "ReferenceIntervalDataProvider")
    public void testQueryAndPrefetch(final Path testReference, final SimpleInterval interval, final String expectedBases ) {
        try (ReferenceDataSource reference = new ReferenceFileSource(testReference))  {
            ReferenceSequence queryResult = reference.queryAndPrefetch(interval);

            Assert.assertEquals(new String(queryResult.getBases()), expectedBases,
                    "Wrong bases returned from queryAndPrefetch() for interval " + interval);
        }
    }

    @Test(dataProvider = "ReferenceIntervalDataProvider")
    public void testQueryAndIterate(final Path testReference, final SimpleInterval interval, final String expectedBases ) {
        try (ReferenceDataSource reference = new ReferenceFileSource(testReference)) {
            Iterator<Byte> queryResultIterator = reference.query(interval);
            List<Byte> queryResult = new ArrayList<>();

            while (queryResultIterator.hasNext()) {
                queryResult.add(queryResultIterator.next());
            }

            Assert.assertEquals(queryResult.size(), expectedBases.length(), "Wrong number of bases returned in iterator from query() call");

            byte[] expectedBytes = expectedBases.getBytes();
            for (int baseIndex = 0; baseIndex < queryResult.size(); ++baseIndex) {
                Assert.assertEquals(queryResult.get(baseIndex).byteValue(), expectedBytes[baseIndex],
                        "Base number " + (baseIndex + 1) + " in iterator from query() call is incorrect");
            }
        }
    }

    /**
     * Test that we can successfully load and query our full-sized B37 reference.
     */
    @Test
    public void testLoadAndQueryB37Reference() {
        try (final ReferenceDataSource ref = new ReferenceFileSource(IOUtils.getPath(b37Reference))) {
            Assert.assertEquals(ref.getSequenceDictionary().getSequences().size(), 85, "Wrong number of contigs in reference sequence dictionary");

            Assert.assertTrue(Arrays.equals(ref.queryAndPrefetch("1", 10000000, 10000005).getBases(), new byte[]{ 'A', 'A', 'C', 'C', 'C', 'C' }), "Wrong reference bases returned for query on 1:10000000-10000005");
        }
    }

    /**
     * Test that we can successfully load and query our full-sized HG38 reference.
     */
    @Test
    public void testLoadAndQueryHG38Reference() {
        try (final ReferenceDataSource ref = new ReferenceFileSource(IOUtils.getPath(hg38Reference))) {
            Assert.assertEquals(ref.getSequenceDictionary().getSequences().size(), 3366, "Wrong number of contigs in reference sequence dictionary");

            Assert.assertTrue(Arrays.equals(ref.queryAndPrefetch("chr1", 10000000, 10000005).getBases(), new byte[]{ 'C', 'A', 'G', 'G', 'T', 'G' }), "Wrong reference bases returned for query on chr1:10000000-10000005");
        }
    }
}
