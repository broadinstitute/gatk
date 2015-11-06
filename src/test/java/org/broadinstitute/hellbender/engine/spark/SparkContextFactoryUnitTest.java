package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.SparkConf;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Map;

public class SparkContextFactoryUnitTest extends BaseTest {

    @Test
    public void testSetupSparkConf(){
        final String appName = "appName";
        final String master = "master";
        final SparkConf sparkConf = SparkContextFactory.setupSparkConf(appName, master, SparkContextFactory.TEST_ATTRIBUTES);
        Assert.assertEquals(sparkConf.get("spark.master"), master);
        Assert.assertEquals(sparkConf.get("spark.app.name"), appName);
        for(Map.Entry<String,String> entry: SparkContextFactory.TEST_ATTRIBUTES.entrySet()){
            Assert.assertEquals(sparkConf.get(entry.getKey()), entry.getValue());
        }
    }

}
