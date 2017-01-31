package org.broadinstitute.hellbender.conf;

import org.apache.commons.configuration.AbstractConfiguration;
import org.apache.commons.configuration.BaseConfiguration;
import org.apache.commons.configuration.Configuration;
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
        return setGATKPropertyList(GATKConf.GATK_PROPERTY_PREFIX, new ArrayList<>(DEFAULT_TOOLS_PACKAGES))
                .setGATKPropertyList(GATKConf.TOOL_CLASS_LIST_KEY, new ArrayList<>())
                .setGATKPropertyList(GATKConf.CODEC_PACKAGE_LIST_KEY, new ArrayList<>(DEFAULT_CODEC_PACKAGES));
    }

    public GATKConf getConfiguration() {
        return new GATKConf(configuration);
    }

    public GATKConfBuilder setGATKPropertyList(final String key, final List<String> value) {
        configuration.setProperty(GATKConf.GATK_PROPERTY_PREFIX + key, value);
        return this;
    }

    public GATKConfBuilder addGATKProperty(final String key, final String value) {
        configuration.addProperty(GATKConf.GATK_PROPERTY_PREFIX + key, value);
        return this;
    }
}
