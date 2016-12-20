package org.broadinstitute.hellbender.conf;

import org.apache.commons.configuration.AbstractConfiguration;
import org.apache.commons.configuration.BaseConfiguration;
import org.apache.commons.configuration.CombinedConfiguration;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.DefaultConfigurationBuilder;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class GATKConfBuilder {

    public static final List<String> DEFAULT_TOOLS_PACKAGES =
            Collections.singletonList("org.broadinstitute.hellbender");

    public static final List<String> DEFAULT_CODEC_PACKAGES = Arrays.asList("htsjdk.variant",
            "htsjdk.tribble",
            "org.broadinstitute.hellbender.utils.codecs");

    private AbstractConfiguration configuration;

    private boolean loadFile = false;

    public GATKConfBuilder() {
        this.configuration = new BaseConfiguration();
    }

    public GATKConfBuilder(final File file) {
        this.configuration = new DefaultConfigurationBuilder(file);
    }

    public static final GATKConf getDefaultConfiguration() {
        return new GATKConfBuilder().setDefaults().getConfiguration();
    }

    public static final GATKConf useConfiguration(final Configuration configuration) {
        return new GATKConf(configuration);
    }

    public GATKConfBuilder setDefaults() {
        return setToolsPackages(new ArrayList<>(DEFAULT_TOOLS_PACKAGES))
                .setClpClasses(new ArrayList<>())
                .setCodecPackages(new ArrayList<>(DEFAULT_CODEC_PACKAGES));
    }

    public GATKConf getConfiguration() {
        return new GATKConf(configuration);
    }

    private GATKConfBuilder setProperty(final String key, final Object value) {
        configuration.setProperty(GATKConf.GATK_PROPERTY_PREFIX + key, value);
        return this;
    }

    private GATKConfBuilder addProperty(final String key, final Object value) {
        configuration.addProperty(GATKConf.GATK_PROPERTY_PREFIX + key, value);
        return this;
    }

    public GATKConfBuilder setToolsPackages(final List<String> toolsPackages) {
        return setProperty(GATKConf.TOOLS_PACKAGES_KEY, toolsPackages);
    }

    public GATKConfBuilder addToolsPackage(final String toolsPackage) {
        return addProperty(GATKConf.TOOLS_PACKAGES_KEY, toolsPackage);
    }

    public GATKConfBuilder setClpClasses(final List<Class<? extends CommandLineProgram>> toolClpClasses) {
        return setProperty(GATKConf.TOOL_CLASSES_KEY, toolClpClasses);
    }

    public GATKConfBuilder addClpClasses(final String toolClpClass) {
        return setProperty(GATKConf.TOOL_CLASSES_KEY, toolClpClass);
    }

    public GATKConfBuilder setCodecPackages(final List<String> codecPackages) {
        return setProperty(GATKConf.CODEC_PACKAGES_KEY, codecPackages);
    }

    public GATKConfBuilder addCodecPackage(final String codecPackage) {
        return addProperty(GATKConf.CODEC_PACKAGES_KEY, codecPackage);
    }
}
