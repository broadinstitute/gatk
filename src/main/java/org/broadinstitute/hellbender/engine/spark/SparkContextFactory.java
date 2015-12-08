package org.broadinstitute.hellbender.engine.spark;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.serializer.KryoSerializer;

import java.util.Map;

/**
 * Manages creation of the Spark context. In particular, for tests a shared global context is used, since Spark does not
 * support multiple concurrent contexts (see https://issues.apache.org/jira/browse/SPARK-2243), and is susceptible to
 * transient errors if contexts are created and stopped in rapid succession.
 */
public final class SparkContextFactory {

    public static final String DEFAULT_SPARK_MASTER = "local[*]";
    private static final boolean SPARK_DEBUG_ENABLED = Boolean.getBoolean("gatk.spark.debug");

    public static final Map<String, String> TEST_ATTRIBUTES = ImmutableMap.<String, String>builder()
            .put("spark.ui.enabled", Boolean.toString(SPARK_DEBUG_ENABLED))
            .put("spark.kryoserializer.buffer.max", "256m")
            .build();

    public static final Map<String, String> MANDATORY_ATTRIBUTES = ImmutableMap.<String,String>builder()
            .put("spark.serializer", KryoSerializer.class.getCanonicalName())
            .put("spark.kryo.registrator", "org.broadinstitute.hellbender.engine.spark.GATKRegistrator")
            .build();

    public static final Map<String, String> CMDLINE_ATTRIBUTES = ImmutableMap.<String, String>builder()
            .put("spark.kryoserializer.buffer.max", "512m")
            .put("spark.driver.maxResultSize", "0")
            .put("spark.driver.userClassPathFirst", "true")
            .put("spark.executor.userClassPathFirst", "true")
            .put("spark.io.compression.codec", "lzf")
            .put("spark.yarn.executor.memoryOverhead", "600")
            .build();

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
    public static synchronized JavaSparkContext getSparkContext(final String appName, final String master) {
        if (testContextEnabled) {
            final JavaSparkContext context = getTestSparkContext();
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
    public static synchronized void stopSparkContext(final JavaSparkContext context) {
        // only call stop for a non-test context
        if (context != testContext) {
            context.stop();
        }
    }


    /**
     * setup a spark context with the given name, master, and settings
     */
    @VisibleForTesting
    static SparkConf setupSparkConf(final String appName, final String master, final Map<String, String> sparkAttributes) {
        final SparkConf sparkConf = new SparkConf().setAppName(appName).setMaster(master);

        MANDATORY_ATTRIBUTES.forEach(sparkConf::set);
        sparkAttributes.forEach(sparkConf::setIfMissing);
        return sparkConf;
    }
    
    private static JavaSparkContext createSparkContext(final String appName, final String master) {
        final SparkConf sparkConf = setupSparkConf(appName, master, CMDLINE_ATTRIBUTES);
        return new JavaSparkContext(sparkConf);
    }

    private static JavaSparkContext createTestSparkContext() {
        final SparkConf sparkConf = setupSparkConf("TestContext", DEFAULT_SPARK_MASTER, TEST_ATTRIBUTES);
        return new JavaSparkContext(sparkConf);
    }
}
