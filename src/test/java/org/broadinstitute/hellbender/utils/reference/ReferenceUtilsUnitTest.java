package org.broadinstitute.hellbender.utils.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.io.*;
import java.util.Arrays;


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

    @DataProvider(name = "testCalculateMD5Data")
    public Object[][] testCalculateMD5Data() {
        return new Object[][]{
                // reference
                new Object[]{ new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta")},
                new Object[]{ new GATKPath("src/test/resources/large/Homo_sapiens_assembly38.20.21.fasta")},
                //new Object[]{ new GATKPath(hg38Reference) },
        };
    }

    @Test(dataProvider = "testCalculateMD5Data")
    public void testCalculateMD5(GATKPath reference){
        try(ReferenceDataSource source = ReferenceDataSource.of(reference.toPath())) {
            final SAMSequenceDictionary dictionary = source.getSequenceDictionary();
            for (SAMSequenceRecord record : dictionary.getSequences()) {
                SimpleInterval interval = new SimpleInterval(record.getSequenceName(), 1, record.getSequenceLength());
                String md5FromDict = record.getMd5();
                String md5Calculated = ReferenceUtils.calculateMD5(reference, interval);
                Assert.assertEquals(md5Calculated, md5FromDict, record.getSequenceName());
            }
        }
    }
}
