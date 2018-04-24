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
                { new SimpleInterval("1", 1, 3), "NNN" },
                { new SimpleInterval("1", 11041, 11045), "GCAAA" },
                { new SimpleInterval("1", 11210, 11220), "CGGTGCTGTGC" },
                { new SimpleInterval("2", 9995, 10005), "NNNNNNCGTAT" },
                { new SimpleInterval("2", 10001, 10080), "CGTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCAC" },
                { new SimpleInterval("2", 10005, 10084), "TCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCACACCC" },
                { new SimpleInterval("2", 15995, 16000), "TGTCAG" }
        };
    }


    @Test(dataProvider = "ReferenceIntervalDataProvider")
    public void testQueryAndPrefetch( final SimpleInterval interval, final String expectedBases ) {
        try (ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE))  {
            ReferenceSequence queryResult = reference.queryAndPrefetch(interval);

            Assert.assertEquals(new String(queryResult.getBases()), expectedBases,
                    "Wrong bases returned from queryAndPrefetch() for interval " + interval);
        }
    }

    @Test(dataProvider = "ReferenceIntervalDataProvider")
    public void testQueryAndIterate( final SimpleInterval interval, final String expectedBases ) {
        try (ReferenceDataSource reference = new ReferenceFileSource(TEST_REFERENCE)) {
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
}
