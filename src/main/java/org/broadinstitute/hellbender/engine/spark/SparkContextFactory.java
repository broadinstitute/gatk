package org.broadinstitute.hellbender.engine.spark;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.serializer.KryoSerializer;

import java.util.Collections;
import java.util.Map;

/**
 * Manages creation of the Spark context. In particular, for tests a shared global context is used, since Spark does not
 * support multiple concurrent contexts (see https://issues.apache.org/jira/browse/SPARK-2243), and is susceptible to
 * transient errors if contexts are created and stopped in rapid succession.
 */
public final class SparkContextFactory {

    public static final String DEFAULT_SPARK_MASTER = "local[*]";
    private static final boolean SPARK_DEBUG_ENABLED = Boolean.getBoolean("gatk.spark.debug");

    /**
     * GATK will not run without these properties
     * They will always be set unless explicitly overridden with {@link org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions#SPARK_PROPERTY_NAME}
     */
    public static final Map<String, String> MANDATORY_PROPERTIES = ImmutableMap.<String,String>builder()
            .put("spark.serializer", KryoSerializer.class.getCanonicalName())
            .put("spark.kryo.registrator", "org.broadinstitute.hellbender.engine.spark.GATKRegistrator")
            .build();

    /**
     * Default properties, GATK may not run if these are set to bad values.
     * These will be set if there were not already set in the environment.
     */
    public static final Map<String, String> DEFAULT_PROPERTIES = ImmutableMap.<String, String>builder()
            .put("spark.kryoserializer.buffer.max", "512m")
            .put("spark.driver.maxResultSize", "0")
            .put("spark.driver.userClassPathFirst", "true")
            .put("spark.executor.userClassPathFirst", "true")
            .put("spark.io.compression.codec", "lzf")
            .put("spark.yarn.executor.memoryOverhead", "600")
            .build();

    public static final Map<String, String> DEFAULT_TEST_PROPERTIES = ImmutableMap.<String, String>builder()
            .put("spark.ui.enabled", Boolean.toString(SPARK_DEBUG_ENABLED))
            .put("spark.kryoserializer.buffer.max", "256m")
            .build();

    private static boolean testContextEnabled;
    private static JavaSparkContext testContext;

    private SparkContextFactory() {}

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
     * @param overridingProperties properties to set on the spark context, overriding any existing value for the same property
     * @param master the Spark master URL
     */
    public static synchronized JavaSparkContext getSparkContext(final String appName, final Map<String, String> overridingProperties, final String master) {
        if (testContextEnabled) {
            final JavaSparkContext context = getTestSparkContext(overridingProperties);
            if (!master.equals(context.master())) {
                throw new IllegalArgumentException(String.format("Cannot reuse spark context " +
                                "with different spark master URL. Existing: %s, requested: %s.",
                        context.master(), master));
            }
            return context;
        }
        return createSparkContext(appName, overridingProperties, master);
    }




    /**
     * Get the test {@link JavaSparkContext} if it has been registered, otherwise returns null.
     */
    public static synchronized JavaSparkContext getTestSparkContext() {
        return getTestSparkContext(Collections.emptyMap());
    }

    /**
     * Get the test {@link JavaSparkContext} if it has been registered, otherwise returns null.
     *
     * @param overridingProperties properties to set on the spark context, possibly overriding values already set
     */
    public static synchronized JavaSparkContext getTestSparkContext(Map<String, String> overridingProperties) {
        if (testContextEnabled && testContext == null) {
            testContext = createTestSparkContext(overridingProperties);
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
     *
     * @param appName human readable name
     * @param master spark master to use
     * @param suggestedProperties properties to set if no values are set for them already
     * @param overridingProperties properties to force to the given value ignoring values already set
     */
    @VisibleForTesting
    static SparkConf setupSparkConf(final String appName, final String master, final Map<String, String> suggestedProperties, final Map<String,String> overridingProperties) {
        final SparkConf sparkConf = new SparkConf().setAppName(appName).setMaster(master);

        suggestedProperties.forEach(sparkConf::setIfMissing);
        MANDATORY_PROPERTIES.forEach(sparkConf::set);
        overridingProperties.forEach(sparkConf::set);

        return sparkConf;
    }
    
    private static JavaSparkContext createSparkContext(final String appName, Map<String,String> overridingProperties, final String master) {
        final SparkConf sparkConf = setupSparkConf(appName, master, DEFAULT_PROPERTIES, overridingProperties);

        return new JavaSparkContext(sparkConf);
    }

    private static JavaSparkContext createTestSparkContext(Map<String, String> overridingProperties) {
        final SparkConf sparkConf = setupSparkConf("TestContext", DEFAULT_SPARK_MASTER, DEFAULT_TEST_PROPERTIES, overridingProperties);
        return new JavaSparkContext(sparkConf);
    }
}
