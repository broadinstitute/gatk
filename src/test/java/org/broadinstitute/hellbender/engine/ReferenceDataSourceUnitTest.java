package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class ReferenceDataSourceUnitTest extends BaseTest {

    private static final File TEST_REFERENCE = new File(hg19MiniReference);

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentReference() {
        ReferenceDataSource refDataSource = new ReferenceDataSource(new File("/foo/bar/nonexistent.fasta"));
    }

    @Test(expectedExceptions = UserException.MissingReferenceFaiFile.class)
    public void testReferenceWithMissingFaiFile() {
        ReferenceDataSource refDataSource = new ReferenceDataSource(new File(publicTestDir + "fastaWithoutFai.fasta"));
    }

    @Test(expectedExceptions = UserException.MissingReferenceDictFile.class)
    public void testReferenceWithMissingDictFile() {
        ReferenceDataSource refDataSource = new ReferenceDataSource(new File(publicTestDir + "fastaWithoutDict.fasta"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullReference() {
        ReferenceDataSource refDataSource = new ReferenceDataSource(null);
    }

    @Test
    public void testGetSequenceDictionary() {
        ReferenceDataSource refDataSource = new ReferenceDataSource(TEST_REFERENCE);
        SAMSequenceDictionary sequenceDictionary = refDataSource.getSequenceDictionary();
        refDataSource.close();

        Assert.assertEquals(sequenceDictionary.size(), 4, "Wrong number of sequences in sequence dictionary returned from refDataSource.getSequenceDictionary()");
        for ( String contig : Arrays.asList("1", "2", "3", "4") ) {
            Assert.assertNotNull(sequenceDictionary.getSequence(contig), "Sequence dictionary returned from refDataSource.getSequenceDictionary() lacks expected contig " + contig);
        }
    }

    @DataProvider(name = "ReferenceIntervalDataProvider")
    public Object[][] getReferenceIntervals() {
        ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE);
        GenomeLocParser genomeLocParser = new GenomeLocParser(reference.getSequenceDictionary());
        reference.close();

        return new Object[][] {
                { genomeLocParser.createGenomeLoc("1", 1, 3), "NNN" },
                { genomeLocParser.createGenomeLoc("1", 11041, 11045), "GCAAA" },
                { genomeLocParser.createGenomeLoc("1", 11210, 11220), "CGGTGCTGTGC" },
                { genomeLocParser.createGenomeLoc("2", 9995, 10005), "NNNNNNCGTAT" },
                { genomeLocParser.createGenomeLoc("2", 10001, 10080), "CGTATCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCAC" },
                { genomeLocParser.createGenomeLoc("2", 10005, 10084), "TCCCACACACCACACCCACACACCACACCCACACACACCCACACCCACACCCACACACACCACACCCACACACCACACCC" },
                { genomeLocParser.createGenomeLoc("2", 15995, 16000), "TGTCAG" }
        };
    }


    @Test(dataProvider = "ReferenceIntervalDataProvider")
    public void testQueryAndPrefetch( final GenomeLoc interval, final String expectedBases ) {
        ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE);
        ReferenceSequence queryResult = reference.queryAndPrefetch(interval);

        Assert.assertEquals(new String(queryResult.getBases()), expectedBases,
                            "Wrong bases returned from queryAndPrefetch() for interval " + interval);

        reference.close();
    }

    @Test(dataProvider = "ReferenceIntervalDataProvider")
    public void testQueryAndIterate( final GenomeLoc interval, final String expectedBases ) {
        ReferenceDataSource reference = new ReferenceDataSource(TEST_REFERENCE);
        Iterator<Byte> queryResultIterator = reference.query(interval);
        List<Byte> queryResult = new ArrayList<Byte>();

        while ( queryResultIterator.hasNext() ) {
            queryResult.add(queryResultIterator.next());
        }

        Assert.assertEquals(queryResult.size(), expectedBases.length(), "Wrong number of bases returned in iterator from query() call");

        byte[] expectedBytes = expectedBases.getBytes();
        for ( int baseIndex = 0; baseIndex < queryResult.size(); ++baseIndex ) {
            Assert.assertEquals(queryResult.get(baseIndex).byteValue(), expectedBytes[baseIndex],
                    "Base number " + (baseIndex + 1) + " in iterator from query() call is incorrect");
        }

        reference.close();
    }
}
