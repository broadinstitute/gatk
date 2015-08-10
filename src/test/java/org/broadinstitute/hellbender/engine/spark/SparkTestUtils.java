package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;

public class SparkTestUtils {
    /**
     * getTestContext creates a (Java) SparkContext and uses the in-memory master and (2) workers.
     * It also does the other stuff we expect (Kryo + kryo registrator). No prerequisites for use.
     * @return new JavaSparkContext.
     */
    public static JavaSparkContext getTestContext() {
        SparkConf sparkConf = new SparkConf().setAppName("TestContext")
                .setMaster("local[2]").set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
                .set("spark.kryo.registrator", "org.broadinstitute.hellbender.engine.spark.GATKRegistrator");
        return new JavaSparkContext(sparkConf);
    }
}
