package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.launcher.SparkLauncher;
import org.apache.spark.serializer.KryoSerializer;

import java.io.File;
import java.io.IOException;

/**
 * Manages creation of the Spark context. In particular, for tests a shared global context is used, since Spark does not
 * support multiple concurrent contexts (see https://issues.apache.org/jira/browse/SPARK-2243), and is susceptible to
 * transient errors if contexts are created and stopped in rapid succession.
 */
public final class SparkContextFactory {

    private static boolean testContextEnabled;
    private static JavaSparkContext testContext;

    private SparkContextFactory() {
    }

    /**
     * Register a shared global {@link JavaSparkContext} for this JVM; this should only be used for testing.
     */
    public static synchronized void enableTestSparkContext() {
        testContextEnabled = true;
    }

    /**
     * Get a {@link JavaSparkContext}. If the test context has been set then it will be returned.
     *
     * @param appName the name of the application to run
     * @param master the Spark master URL
     */
    public static synchronized JavaSparkContext getSparkContext(String appName, String master) {
        if (testContextEnabled) {
            JavaSparkContext context = getTestSparkContext();
            if (!master.equals(context.master())) {
                throw new IllegalArgumentException(String.format("Cannot reuse spark context " +
                                "with different spark master URL. Existing: %s, requested: %s.",
                        context.master(), master));
            }
            return context;
        }
        return createSparkContext(appName, master);
    }

    /**
     * Get the test {@link JavaSparkContext} if it has been registered, otherwise returns null.
     *
     */
    public static synchronized JavaSparkContext getTestSparkContext() {
        if (testContextEnabled && testContext == null) {
            testContext = createTestSparkContext();
            Runtime.getRuntime().addShutdownHook(new Thread() {
                @Override
                public void run() {
                    testContext.stop();
                }
            });
        }
        return testContext;
    }

    /**
     * Stop a {@link JavaSparkContext}, unless it is the test context.
     *
     * @param context the context to stop
     */
    public static synchronized void stopSparkContext(JavaSparkContext context) {
        // only call stop for a non-test context
        if (context != testContext) {
            context.stop();
        }
    }

    private static JavaSparkContext createSparkContext(String appName, String master) {
        SparkConf sparkConf = new SparkConf().setAppName(appName)
                .setMaster(master)
                .set("spark.serializer", KryoSerializer.class.getCanonicalName())
                .set("spark.kryo.registrator", "org.broadinstitute.hellbender.engine.spark.GATKRegistrator")
                .set("spark.kryoserializer.buffer.max", "512m");

        return new JavaSparkContext(sparkConf);
    }

    private static JavaSparkContext createTestSparkContext() {
        String path = null;
        String jarLibPath = null;
        try {
            path = "/Users/davidada/apps/spark-1.4.1-bin-hadoop2.6/lib/*:" + new File("build/classes/main").getCanonicalPath();
            jarLibPath = new File("build/lib/*").getCanonicalPath();
        } catch (IOException e) {
            e.printStackTrace();
        }
        SparkConf sparkConf = new SparkConf().setAppName("TestContext")
                .setSparkHome("/Users/davidada/apps/spark-1.4.1-bin-hadoop2.6")
                //.set("spark.driver.userClassPathFirst", "true")
                //.set("spark.executor.userClassPathFirst", "true")
        .set(SparkLauncher.EXECUTOR_EXTRA_CLASSPATH, path + ":" + jarLibPath)
                .set(SparkLauncher.DRIVER_EXTRA_CLASSPATH, path + ":" + jarLibPath)
                .setMaster("spark://davidada-macbookpro2.roam.corp.google.com:7077") //.setMaster("local[2]")
                .set("spark.serializer", KryoSerializer.class.getCanonicalName())
                .set("spark.kryo.registrator", "org.broadinstitute.hellbender.engine.spark.GATKRegistrator")
                ;//.set("spark.ui.enabled", "false");

        return new JavaSparkContext(sparkConf);
    }
}
