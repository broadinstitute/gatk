package org.broadinstitute.hellbender.conf;

import org.apache.commons.configuration.Configuration;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;

/**
 * Configuration for GATK. Properties are prefixed by {@link #GATK_PROPERTY_PREFIX}.
 * <p>
 * - {@link #TOOL_PACKAGE_LIST_KEY}: list of java packages in which to search for classes that extend CommandLineProgram.
 * - {@link #TOOL_CLASS_LIST_KEY}: list of single classes to include (e.g. required input pre-processing tools).
 * - {@link #CODEC_PACKAGE_LIST_KEY}: packages to search at startup to look for FeatureCodecs.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class GATKConf {

    public static final String GATK_PROPERTY_PREFIX = "gatk.";

    /**
     * The packages we wish to include in our command line.
     */
    public final static String TOOL_PACKAGE_LIST_KEY = "toolsPackages";

    /**
     * The single classes we wish to include in our command line.
     */
    public final static String TOOL_CLASS_LIST_KEY = "toolClasses";

    /**
     * We will search these packages at startup to look for FeatureCodecs.
     */
    public final static String CODEC_PACKAGE_LIST_KEY = "codecPackages";

    final Configuration configuration;

    GATKConf(final Configuration configuration) {
        Utils.nonNull(configuration);
        this.configuration = configuration;
    }

    @SuppressWarnings("unchecked")
    public List<String> getGATKProperty(final String propertyKey) {
        return Collections.unmodifiableList(configuration.getList(GATK_PROPERTY_PREFIX + propertyKey));
    }
}
