package org.broadinstitute.hellbender.conf;

import org.apache.commons.configuration.Configuration;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Configuration for GATK. Properties are prefixed by {@link #GATK_PROPERTY_PREFIX}.
 *
 * - {@link #TOOLS_PACKAGES_KEY}: list of java packages in which to search for classes that extend CommandLineProgram.
 * - {@link #TOOL_CLASSES_KEY}: list of single classes to include (e.g. required input pre-processing tools).
 * - {@link #CODEC_PACKAGES_KEY}: packages to search at startup to look for FeatureCodecs.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class GATKConf {

    public static final String GATK_PROPERTY_PREFIX = "gatk.";

    public final static String TOOLS_PACKAGES_KEY = "toolsPackages";
    public final static String TOOL_CLASSES_KEY = "toolClasses";
    public final static String CODEC_PACKAGES_KEY = "codecPackages";

    final Configuration configuration;

    GATKConf(final Configuration configuration) {
        Utils.nonNull(configuration);
        this.configuration = configuration;
    }

    /**
     * The packages we wish to include in our command line.
     */
    @SuppressWarnings("unchecked")
    public List<String> getToolPackages() {
        return configuration.getList(GATK_PROPERTY_PREFIX +  TOOLS_PACKAGES_KEY);
    }

    /**
     * The single classes we wish to include in our command line.
     */
    @SuppressWarnings("unchecked")
    public List<Class<? extends CommandLineProgram>> getToolClasses() {
        return configuration.getList(GATK_PROPERTY_PREFIX +  TOOL_CLASSES_KEY);
    }

    /**
     * We will search these packages at startup to look for FeatureCodecs.
     */
    @SuppressWarnings("unchecked")
    public List<String> getCodecPackages() {
        return configuration.getList(GATK_PROPERTY_PREFIX +  CODEC_PACKAGES_KEY);
    }

}
