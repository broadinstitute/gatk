package org.broadinstitute.hellbender.utils.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.nio.file.Paths;
import java.util.*;

import java.io.*;


public class ReferenceUtilsUnitTest extends GATKBaseTest {

    @Test
    public void testLoadFastaDictionaryFromFile() {
        final File referenceDictionaryFile = new File(ReferenceUtils.getFastaDictionaryFileName(hg19MiniReference));
        final SAMSequenceDictionary dictionary = ReferenceUtils.loadFastaDictionary(referenceDictionaryFile);

        Assert.assertNotNull(dictionary, "Sequence dictionary null after loading");
        Assert.assertEquals(dictionary.size(), 4, "Wrong sequence dictionary size after loading");
    }

    private static final class ClosingAwareFileInputStream extends FileInputStream {
        private boolean isClosed;

        public ClosingAwareFileInputStream( final File file ) throws FileNotFoundException {
            super(file);
            isClosed = false;
        }

        @Override
        public void close() throws IOException {
            super.close();
            isClosed = true;
        }

        public boolean isClosed() {
            return isClosed;
        }
    }

    @Test
    public void testLoadFastaDictionaryFromStream() throws IOException {
        try ( final ClosingAwareFileInputStream referenceDictionaryStream = new ClosingAwareFileInputStream(new File(ReferenceUtils.getFastaDictionaryFileName(hg19MiniReference))) ) {
            final SAMSequenceDictionary dictionary = ReferenceUtils.loadFastaDictionary(referenceDictionaryStream);

            Assert.assertNotNull(dictionary, "Sequence dictionary null after loading");
            Assert.assertEquals(dictionary.size(), 4, "Wrong sequence dictionary size after loading");

            Assert.assertFalse(referenceDictionaryStream.isClosed(), "InputStream was improperly closed by ReferenceUtils.loadFastaDictionary()");
        }
    }

    @DataProvider
    public Object[][] testData_testLoadWrongFormatFastaDictionary() {

        // Trying to ingest a fasta file as a dict file to induce a MalformedFile Exception.
        return new Object[][] {
                {hg19MiniReference, true},
                {hg19MiniReference, false}
        };
    }

    @Test(expectedExceptions = UserException.MalformedFile.class, dataProvider = "testData_testLoadWrongFormatFastaDictionary")
    public void testLoadWrongFormatFastaDictionary(final String fastaFileName, boolean loadAsStream ) throws IOException {

        File testFastaFile = new File(fastaFileName);

        if ( loadAsStream ) {
            // Test loadFastaDictionary with Stream argument:
            try ( final ClosingAwareFileInputStream referenceDictionaryStream = new ClosingAwareFileInputStream(testFastaFile) ) {
                ReferenceUtils.loadFastaDictionary(referenceDictionaryStream);
            }
        }
        else {
            // Test loadFastaDictionary with File argument:
            ReferenceUtils.loadFastaDictionary(testFastaFile);
        }
    }

    @Test(groups = {"bucket"})
    public void testLoadFastaDictionaryFromGCSBucket() throws IOException {
        final String bucketDictionary = getGCPTestInputPath() + "org/broadinstitute/hellbender/utils/ReferenceUtilsTest.dict";

        try ( final InputStream referenceDictionaryStream = BucketUtils.openFile(bucketDictionary) ) {
            final SAMSequenceDictionary dictionary = ReferenceUtils.loadFastaDictionary(referenceDictionaryStream);

            Assert.assertNotNull(dictionary, "Sequence dictionary null after loading");
            Assert.assertEquals(dictionary.size(), 4, "Wrong sequence dictionary size after loading");
        }
    }

    @Test(dataProvider =  "dumpReferenceSourceData")
    public void testDumpReferenceSourceInFiles(final ReferenceSource source, final int basesPerLine)
        throws IOException
    {
        final File whereTo = createTempFile("test-dump-ref-source-", ".fasta");
        final File whereToIndex = new File(whereTo.getParent(), whereTo.getName() + ".fai");
        whereToIndex.deleteOnExit();
        final ReferenceSource resultSource = ReferenceUtils.dumpReferenceSource(source, Paths.get(whereTo.toURI()), basesPerLine);
        Assert.assertNotNull(resultSource);
        Assert.assertTrue(whereTo.isFile());
        // Check on the index:
        Assert.assertEquals(resultSource.getReferenceSequenceDictionary(), source.getReferenceSequenceDictionary());
    }


    @DataProvider(name = "dumpReferenceSourceData")
    public Object[][] dumpReferenceSourceData()
        throws IOException
    {
        final RandomDNA randomDNA = new RandomDNA();
        final File TWO_CHR_10k_4k5;
        randomDNA.nextFasta(TWO_CHR_10k_4k5 = createTempFile("test-ref", ".fasta"),
                new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("seq1", 10_000),
                        new SAMSequenceRecord("seq2", 4_500))), 76);
        final FastaSequenceIndex index = FastaSequenceIndexCreator.buildFromFasta(Paths.get(TWO_CHR_10k_4k5.toURI()));
        final File TWO_CHR_10k_4k5_index = new File(TWO_CHR_10k_4k5.getParent(), TWO_CHR_10k_4k5.getName() + ".fai");
        index.write(TWO_CHR_10k_4k5_index.toPath());
        TWO_CHR_10k_4k5_index.deleteOnExit();
        return new Object[][] {
                new Object[] { new ReferenceFileSource(TWO_CHR_10k_4k5.toPath()) , 56 }
        };
    }
}
