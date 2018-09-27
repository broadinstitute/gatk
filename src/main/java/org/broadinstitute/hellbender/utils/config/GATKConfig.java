package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Accessible;
import org.aeonbits.owner.Mutable;
import org.aeonbits.owner.Config.LoadPolicy;
import org.aeonbits.owner.Config.LoadType;
import org.aeonbits.owner.Config.Sources;

import java.util.List;

/**
 * Configuration file for GATK options.
 * All specified {@code Sources} will be loaded.
 * The {@link LoadPolicy} is set to {@link LoadType#MERGE}, which specifies that if a configuration option is not found
 * in the first source, the option will be sought in all following sources until a definition is found.
 * If the option is not specified in any file, the coded default value will be used (as defined by @DefaultValue).
 *
 * The load order is always resolved "top-down" by declaration order in the @Sources annotation.
 *
 * In this case, the load order is:
 *        1)   "file:${" + GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME + "}",
 *        2)   "classpath:${" + GATKConfig.CONFIG_FILE_VARIABLE_CLASS_PATH + "}",
 *        3)   "file:GATKConfig.properties",
 *        4)   "classpath:org/broadinstitute/hellbender/utils/config/GATKConfig.properties"
 *        5)   hard-coded values specified by @DefaultValue
 *
 */
@LoadPolicy(LoadType.MERGE)
@Sources({
        "file:${" + GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME + "}",                 // Variable for file loading
        "classpath:${" + GATKConfig.CONFIG_FILE_VARIABLE_CLASS_PATH + "}",           // Variable for class path loading
        "file:GATKConfig.properties",                                                // Default path
        "classpath:org/broadinstitute/hellbender/utils/config/GATKConfig.properties" // Class path
})
public interface GATKConfig extends Mutable, Accessible {

    // =================================================================================================================
    // =================================================================================================================
    // Meta Options:
    // =================================================================================================================
    // =================================================================================================================

    /**
     * Name of the configuration file variable to be used in the {@link Sources} annotation for {@link GATKConfig}
     * as a place to find the configuration file corresponding to this interface.
     */
    String CONFIG_FILE_VARIABLE_FILE_NAME = "GATKConfig.pathToGatkConfig";

    /**
     * Name of the configuration file variable to be used in the {@link Sources} annotation for {@link GATKConfig}
     * as a place to find the configuration file corresponding to this interface.
     */
    String CONFIG_FILE_VARIABLE_CLASS_PATH = "GATKConfig.classPathToGatkConfig";

    // =================================================================================================================
    // =================================================================================================================
    // System Options:
    // =================================================================================================================
    // =================================================================================================================

    // ----------------------------------------------------------
    // Miscellaneous Options:
    // ----------------------------------------------------------

    @SystemProperty
    @Key("gatk_stacktrace_on_user_exception")
    @DefaultValue("false")
    boolean gatk_stacktrace_on_user_exception();

    // ----------------------------------------------------------
    // HTSJDK Options:
    // ----------------------------------------------------------

    @SystemProperty
    @Key("samjdk.use_async_io_read_samtools")
    @ConverterClass(CustomBooleanConverter.class)
    @DefaultValue("false")
    Boolean samjdk_use_async_io_read_samtools();

    @SystemProperty
    @Key("samjdk.use_async_io_write_samtools")
    @DefaultValue("true")
    boolean samjdk_use_async_io_write_samtools();

    @SystemProperty
    @Key("samjdk.use_async_io_write_tribble")
    @DefaultValue("false")
    boolean samjdk_use_async_io_write_tribble();

    @SystemProperty
    @Key("samjdk.compression_level")
    @DefaultValue("2")
    int samjdk_compression_level();

    // ----------------------------------------------------------
    // Spark Options:
    // ----------------------------------------------------------

    @SystemProperty
    @Key("spark.kryoserializer.buffer.max")
    @DefaultValue("512m")
    String spark_kryoserializer_buffer_max();

    @SystemProperty
    @Key("spark.driver.maxResultSize")
    @DefaultValue("0")
    int spark_driver_maxResultSize();

    @SystemProperty
    @Key("spark.driver.userClassPathFirst")
    @DefaultValue("true")
    boolean spark_driver_userClassPathFirst();

    @SystemProperty
    @Key("spark.io.compression.codec")
    @DefaultValue("lzf")
    String spark_io_compression_codec();

    @SystemProperty
    @Key("spark.yarn.executor.memoryOverhead")
    @DefaultValue("600")
    int spark_yarn_executor_memoryOverhead();

    @SystemProperty
    @Key("spark.driver.extraJavaOptions")
    @DefaultValue("")
    String spark_driver_extraJavaOptions();

    @SystemProperty
    @Key("spark.executor.extraJavaOptions")
    @DefaultValue("")
    String spark_executor_extraJavaOptions();

    // =================================================================================================================
    // =================================================================================================================
    // GATK  Options:
    // =================================================================================================================
    // =================================================================================================================
    
    // ----------------------------------------------------------
    // Miscellaneous Options:
    // ----------------------------------------------------------

    @DefaultValue("htsjdk.variant,htsjdk.tribble,org.broadinstitute.hellbender.utils.codecs")
    List<String> codec_packages();

    // ----------------------------------------------------------
    // GATKTool Options:
    // ----------------------------------------------------------

    @DefaultValue("40")
    int cloudPrefetchBuffer();

    @DefaultValue("-1")
    int cloudIndexPrefetchBuffer();

    @DefaultValue("20")
    int gcsMaxRetries();

    @DefaultValue("")
    String gcsProjectForRequesterPays();

    @DefaultValue("true")
    boolean createOutputBamIndex();
}
