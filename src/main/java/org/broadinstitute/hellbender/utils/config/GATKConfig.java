package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Accessible;
import org.aeonbits.owner.Mutable;
import org.aeonbits.owner.Config.LoadPolicy;
import org.aeonbits.owner.Config.LoadType;
import org.aeonbits.owner.Config.Sources;

import java.util.List;

/**
 * Configuration file for GATK options.
 */
@LoadPolicy(LoadType.MERGE)
@Sources({
        "file:${" + GATKConfig.CONFIG_FILE_VARIABLE_NAME + "}",                         // Variable for file loading
        "file:GATKConfig.properties",                                                   // Default path
        "classpath:org/broadinstitute/hellbender/utils/config/GATKConfig.properties"    // Class path
})
public interface GATKConfig extends Mutable, Accessible {

    // =================================================================================================================
    // =================================================================================================================
    // Meta Options:
    // =================================================================================================================
    // =================================================================================================================

    String CONFIG_FILE_VARIABLE_NAME = "pathToGatkConfig";

    // =================================================================================================================
    // =================================================================================================================
    // System Options:
    // =================================================================================================================
    // =================================================================================================================

    String SYSTEM_PROPERTY_PREFIX = "system.";

    // ----------------------------------------------------------
    // General Options:
    // ----------------------------------------------------------

    @Key("system.gatk_stacktrace_on_user_exception")
    @DefaultValue("true")
    boolean gatk_stacktrace_on_user_exception();

    // ----------------------------------------------------------
    // SAMJDK Options:
    // ----------------------------------------------------------

    @Key("system.samjdk.use_async_io_read_samtools")
    @ConverterClass(CustomBooleanConverter.class)
    @DefaultValue("false")
    Boolean samjdk_use_async_io_read_samtools();

    @Key("system.samjdk.use_async_io_write_samtools")
    @DefaultValue("true")
    boolean samjdk_use_async_io_write_samtools();

    @Key("system.samjdk.use_async_io_write_tribble")
    @DefaultValue("false")
    boolean samjdk_use_async_io_write_tribble();

    @Key("system.samjdk.compression_level")
    @DefaultValue("1")
    int samjdk_compression_level();

    // ----------------------------------------------------------
    // Spark Options:
    // ----------------------------------------------------------

    @Key("system.spark.kryoserializer.buffer.max")
    @DefaultValue("512m")
    String spark_kryoserializer_buffer_max();

    @Key("system.spark.driver.maxResultSize")
    @DefaultValue("0")
    int spark_driver_maxResultSize();

    @Key("system.spark.driver.userClassPathFirst")
    @DefaultValue("true")
    boolean spark_driver_userClassPathFirst();

    @Key("system.spark.io.compression.codec")
    @DefaultValue("lzf")
    String spark_io_compression_codec();

    @Key("system.spark.yarn.executor.memoryOverhead")
    @DefaultValue("600")
    int spark_yarn_executor_memoryOverhead();

    @Key("system.spark.driver.extraJavaOptions")
    @DefaultValue("")
    String spark_driver_extraJavaOptions();

    @Key("system.spark.executor.extraJavaOptions")
    @DefaultValue("")
    String spark_executor_extraJavaOptions();

    // ----------------------------------------------------------
    // Other Options:
    // ----------------------------------------------------------

    @Key("system.snappy.disable")
    @DefaultValue("true")
    boolean snappy_disable();

    // =================================================================================================================
    // =================================================================================================================
    // System Options:
    // =================================================================================================================
    // =================================================================================================================
    
    // ----------------------------------------------------------
    // General Options:
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
    int gcsMaxNumRetries();

    @DefaultValue("true")
    boolean createOutputBamIndex();
}