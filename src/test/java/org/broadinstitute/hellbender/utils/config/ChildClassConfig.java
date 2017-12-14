package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Accessible;
import org.aeonbits.owner.Config;
import org.aeonbits.owner.Mutable;

import java.util.List;

/**
 * Simple test configuration file to check overrides and other configuration features.
 */
@Config.LoadPolicy(Config.LoadType.MERGE)
@Config.Sources({
        "file:/etc/jaiow",                          // Test of non-existent bad path
        "classpath:org/broadinstitute/hellbender/utils/config/ChildClassConfig.properties",
})
public interface ChildClassConfig extends BasicTestConfigWithClassPathOverridesAndVariableFile {

    @DefaultValue("true")
    @ConverterClass(CustomBooleanConverter.class)
    Boolean newCustomBooleanThatDefaultsToTrue();
}
