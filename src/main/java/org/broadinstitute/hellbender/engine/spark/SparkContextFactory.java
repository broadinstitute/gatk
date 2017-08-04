package org.broadinstitute.hellbender.engine.spark;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.serializer.KryoSerializer;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Map;

/**
 * Manages creation of the Spark context. In particular, for tests a shared global context is used, since Spark does not
 * support multiple concurrent contexts (see https://issues.apache.org/jira/browse/SPARK-2243), and is susceptible to
 * transient errors if contexts are created and stopped in rapid succession.
 */
public final class SparkContextFactory {

    public static final String DEFAULT_SPARK_MASTER = determineDefaultSparkMaster();
    private static final boolean SPARK_DEBUG_ENABLED = Boolean.getBoolean("gatk.spark.debug");
    private static final String SPARK_CORES_ENV_VARIABLE = "GATK_TEST_SPARK_CORES";
    private static final String TEST_PROJECT_ENV_VARIABLE = "HELLBENDER_TEST_PROJECT";
    private static final String TEST_JSON_KEYFILE_ENV_VARIABLE = "HELLBENDER_JSON_SERVICE_ACCOUNT_KEY";

    private static final Logger logger = LogManager.getLogger(SparkContextFactory.class);

    /**
     * GATK will not run without these properties
     * They will always be set unless explicitly overridden with {@link org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions#SPARK_PROPERTY_NAME}
     */
    public static final Map<String, String> MANDATORY_PROPERTIES = ImmutableMap.<String,String>builder()
            .put("spark.serializer", KryoSerializer.class.getCanonicalName())
            .put("spark.kryo.registrator", GATKRegistrator.class.getCanonicalName())
            // remap the Hadoop FS implementation for file:// URIs to avoid writing CRC files for local files
            // note that we don't use Hadoop's RawLocalFileSystem since it doesn't extend LocalFileSystem
            .put("spark.hadoop.fs.file.impl", NonChecksumLocalFileSystem.class.getCanonicalName())
            .build();

    /**
     * Default properties, GATK may not run if these are set to bad values.
     * These will be set if there were not already set in the environment.
     */
    public static final Map<String, String> DEFAULT_PROPERTIES = ImmutableMap.<String, String>builder()
            .put("spark.kryoserializer.buffer.max", "512m")
            .put("spark.driver.maxResultSize", "0")
            .put("spark.driver.userClassPathFirst", "true")
            .put("spark.io.compression.codec", "lzf")
            .put("spark.yarn.executor.memoryOverhead", "600")
            .build();

    public static final Map<String, String> DEFAULT_TEST_PROPERTIES = ImmutableMap.<String, String>builder()
            .put("spark.ui.enabled", Boolean.toString(SPARK_DEBUG_ENABLED))
            .put("spark.kryoserializer.buffer.max", "256m")
            .put("spark.hadoop.fs.file.impl.disable.cache", "true") // so NonChecksumLocalFileSystem is not cached between tests
            .putAll(getGcsHadoopAdapterTestProperties())
            .build();

    /**
     * @return checks if the necessary environment variables are present in order to configure the gcs-hadoop adapter
     * and returns a map containing the configuration if they are available
     * returns an empty map otherwise
     */
    private static Map<String,String> getGcsHadoopAdapterTestProperties(){
        final String testProject = System.getenv(TEST_PROJECT_ENV_VARIABLE);
        final String testKeyFile = System.getenv(TEST_JSON_KEYFILE_ENV_VARIABLE);
        if( testProject == null || testKeyFile == null) {
            logger.warn("Environment variables " + TEST_PROJECT_ENV_VARIABLE + " and " + TEST_JSON_KEYFILE_ENV_VARIABLE +
                                " must be set or the GCS hadoop connector will not be configured properly");
            return Collections.emptyMap();
        } else {
            return ImmutableMap.<String, String>builder()
                    .put("spark.hadoop.fs.gs.impl", "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem")
                    .put("spark.hadoop.fs.AbstractFileSystem.gs.impl", "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS")
                    .put("spark.hadoop.fs.gs.project.id", testProject)
                    .put("spark.hadoop.google.cloud.auth.service.account.json.keyfile", testKeyFile)
                    .build();
        }
    }

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
            Utils.validateArg(master.equals(context.master()), () -> String.format("Cannot reuse spark context " +
                            "with different spark master URL. Existing: %s, requested: %s.", context.master(), master));
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

    /**
     * Create the default Spark master, determines the number of cores it should use. Applicable to Spark test only.
     *   Read the specification from the environmental variable GATK_TEST_SPARK_CORES
     *      If the enviromental variable is not set,  use all available cores as in "local[*]"
     *      If the value is a positive integer, use the value
     *      If the value is invalid (strings, empty, etc), throw an UserException
     *      If the value is a negative interger or zero, throw an UserException
     */
    private static String determineDefaultSparkMaster() {
	final String defaultSparkMasterString = "local[*]";
	String sparkMasterString;

	String sparkSpecFromEnvironment = System.getenv( SPARK_CORES_ENV_VARIABLE );
	if ( null == sparkSpecFromEnvironment ) {
	    sparkMasterString = defaultSparkMasterString;
	} else {
	    int numSparkCoresFromEnv = 0;
	    try {
		numSparkCoresFromEnv = Integer.parseInt( sparkSpecFromEnvironment );
	    } catch ( NumberFormatException e ) {
		throw new UserException("Illegal number of cores specified in " + SPARK_CORES_ENV_VARIABLE + ". Positive integers only");
	    }
	    
	    if ( numSparkCoresFromEnv > 0 ) {
		sparkMasterString = String.format("local[%d]", numSparkCoresFromEnv );
	    } else {
		throw new UserException("Illegal number of cores specified in " + SPARK_CORES_ENV_VARIABLE + ". Number of cores must be positive");
	    }
	}
	return sparkMasterString;
    }
}
