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
        "file:${" + GATKConfig.CONFIG_FILE_VARIABLE_NAME + "}",
        "file:/etc/jaiow",                          // Test of non-existent bad path
        "classpath:org/broadinstitute/hellbender/utils/config/BasicTestConfigWithClassPathOverrides.properties",
})
public interface BasicTestConfigWithClassPathOverridesAndVariableFile extends Mutable, Accessible {

    @DefaultValue("true")
    boolean booleanDefTrue();

    @DefaultValue("false")
    boolean booleanDefFalse();

    @DefaultValue("207")
    int intDef207();

    @DefaultValue("string1,string2,string3,string4")
    List<String> listOfStringTest();
}
