package org.broadinstitute.hellbender.utils.dataflow;

import org.testng.Assert;
import org.testng.annotations.Test;

public class BucketUtilsTest {

    @Test
    public void testIsCloudStorageURL(){
        Assert.assertTrue(BucketUtils.isCloudStorageUrl("gs://a_bucket/bucket"));
        Assert.assertFalse(BucketUtils.isCloudStorageUrl("localFile"));
    }
}