package org.broadinstitute.hellbender.utils.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

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
}
