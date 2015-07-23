package org.broadinstitute.hellbender.utils.reference;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;


public class ReferenceUtilsUnitTest extends BaseTest {

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

    @Test(groups = {"bucket"})
    public void testLoadFastaDictionaryFromGCSBucket() throws IOException {
        final String bucketDictionary = getDataflowTestInputPath() + "org/broadinstitute/hellbender/utils/ReferenceUtilsTest.dict";
        final GCSOptions popts = PipelineOptionsFactory.as(GCSOptions.class);
        popts.setApiKey(getDataflowTestApiKey());

        try ( final InputStream referenceDictionaryStream = BucketUtils.openFile(bucketDictionary, popts) ) {
            final SAMSequenceDictionary dictionary = ReferenceUtils.loadFastaDictionary(referenceDictionaryStream);

            Assert.assertNotNull(dictionary, "Sequence dictionary null after loading");
            Assert.assertEquals(dictionary.size(), 4, "Wrong sequence dictionary size after loading");
        }
    }
}
