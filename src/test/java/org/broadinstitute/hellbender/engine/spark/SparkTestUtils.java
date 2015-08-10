package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;

public class SparkTestUtils {
    public static JavaSparkContext getTestContext() {
        SparkConf sparkConf = new SparkConf().setAppName("TestContext")
                .setMaster("local[2]").set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
                .set("spark.kryo.registrator", "org.broadinstitute.hellbender.engine.spark.GATKRegistrator");
        return new JavaSparkContext(sparkConf);
    }
}
