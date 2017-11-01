package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.ImmutableMap;
import org.apache.spark.SparkConf;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class SparkContextFactoryUnitTest extends GATKBaseTest {

    private static final String prop1 = "spark.value1";
    private static final String value1 = "10";
    private static final String value1_override = "4";

    private static final String prop2 = "spark.value2";
    private static final String value2 = "something";

    private static final String prop3 = "spark.value3";
    private static final String value3 = "something";

    private static final ImmutableMap<String, String> OVERRIDE = new ImmutableMap.Builder<String,String>()
            .put(prop1, value1_override)
            .put(prop2, value2)
            .build();

    private static final ImmutableMap<String, String> SUGGESTED = new ImmutableMap.Builder<String,String>()
            .put(prop1, value1)
            .put(prop3, value3)
            .build();

    @Test(groups = "spark")
    public void testSetupSparkConf(){
        final String appName = "appName";
        final String master = "master";
        final SparkConf sparkConf = SparkContextFactory.setupSparkConf(appName, master, SUGGESTED, OVERRIDE);
        Assert.assertEquals(sparkConf.get("spark.master"), master);
        Assert.assertEquals(sparkConf.get("spark.app.name"), appName);

        Assert.assertEquals(sparkConf.get(prop1), value1_override);
        Assert.assertEquals(sparkConf.get(prop2), value2);
        Assert.assertEquals(sparkConf.get(prop3), value3);
    }



}
