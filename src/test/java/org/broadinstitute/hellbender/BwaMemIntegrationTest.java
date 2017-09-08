package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.*;

public class BwaMemIntegrationTest extends BaseTest {

    private BwaMemIndex index;

    @BeforeClass
    public void loadIndex() {
        final String imageFile = createTempFile(b37_reference_20_21, ".img").toString();
        BwaMemIndex.createIndexImageFromFastaFile(b37_reference_20_21, imageFile);
        index = new BwaMemIndex(imageFile);
    }

    @AfterClass
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
