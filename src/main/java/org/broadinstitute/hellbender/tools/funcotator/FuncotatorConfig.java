package org.broadinstitute.hellbender.tools.funcotator;

import org.aeonbits.owner.Accessible;
import org.aeonbits.owner.Config.LoadPolicy;
import org.aeonbits.owner.Config.LoadType;
import org.aeonbits.owner.Config.Sources;
import org.aeonbits.owner.Mutable;
import org.broadinstitute.hellbender.utils.config.CustomBooleanConverter;

import java.util.List;
import java.util.Set;

/**
 * Configuration file for Funcotator options.
 * All specified {@code Sources} will be loaded.
 * The {@link LoadPolicy} is set to {@link LoadType#MERGE}, which specifies that if a configuration option is not found
 * in the first source, the option will be sought in all following sources until a definition is found.
 * If the option is not specified in any file, the coded default value will be used (as defined by @DefaultValue).
 *
 * The load order is always resolved "top-down" by declaration order in the @Sources annotation.
 *g
 * In this case, the load order is:
 *        1)   "file:${" + GATKConfig.CONFIG_FILE_VARIABLE_FILE_NAME + "}",
 *
 */
@LoadPolicy(LoadType.MERGE)
@Sources({
        "file:${" + FuncotatorConfig.CONFIG_FILE_VARIABLE_FILE_NAME + "}",                 // Variable for file loading
})
public interface FuncotatorConfig extends Mutable, Accessible {

    // =================================================================================================================
    // =================================================================================================================
    // Meta Options:
    // =================================================================================================================
    // =================================================================================================================

    /**
     * Name of the configuration file variable to be used in the {@link Sources} annotation for {@link FuncotatorConfig}
     * as a place to find the configuration file corresponding to this interface.
     */
    String CONFIG_FILE_VARIABLE_FILE_NAME = "FuncotatorConfig.pathToGatkConfig";

    // =================================================================================================================
    // =================================================================================================================
    // Runtime options:
    // =================================================================================================================
    // =================================================================================================================

    @Key(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME)
    List<String> dataSourceDirectories();

    @Key(FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME)
    @ConverterClass(CustomBooleanConverter.class)
    @DefaultValue("" + FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_DEFAULT_VALUE)
    Boolean removeFilteredVariants();

    // TODO: Find a way to use the default value for this as defined in FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE
    @Key(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME)
    @DefaultValue("CANONICAL")
    TranscriptSelectionMode transcriptSelectionMode();

    @Key(FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME)
    @DefaultValue("")
    Set<String> transcriptList();

    @Key(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME)
    @DefaultValue("")
    List<String> annotationDefaults();

    @Key(FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME)
    @DefaultValue("")
    List<String> annotationOverrides();

    @Key(FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_ARG_NAME)
    @DefaultValue("" + FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_DEFAULT_VALUE)
    int lookaheadCacheSize();
}
