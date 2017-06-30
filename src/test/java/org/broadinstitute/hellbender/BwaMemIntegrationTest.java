package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

public class BwaMemIntegrationTest extends BaseTest {

    private BwaMemIndex index;

    @BeforeSuite
    public void loadIndex() {
        final String imageFile = createTempFile(TestResources.b37_reference_20_21, ".img").toString();
        BwaMemIndex.createIndexImage(TestResources.b37_reference_20_21, imageFile);
        index = new BwaMemIndex(imageFile);
    }

    @AfterSuite
    public void unloadIndex() {
        index.close();
        index = null;
    }

    @Test
    public void testAlignSingleRead() throws Exception {
        BwaMemTestUtils.assertCorrectSingleReadAlignment(index);
    }

    @Test
    public void testAlignChimericContig() throws Exception {
        BwaMemTestUtils.assertCorrectChimericContigAlignment(index);
    }

}
