package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Accessible;
import org.aeonbits.owner.Config;
import org.aeonbits.owner.Mutable;

import java.util.List;

/**
 * Simple test configuration file to check overrides and other configuration features.
 */
@Config.LoadPolicy(Config.LoadType.MERGE)
public interface SystemTestConfig extends Mutable, Accessible {

    String CONFIG_FILE_NAME = "SystemTestConfig.config";

    @DefaultValue("true")
    boolean booleanDefTrue();

    @DefaultValue("false")
    boolean booleanDefFalse();

    @DefaultValue("207")
    int intDef207();

    @DefaultValue("string1,string2,string3,string4")
    List<String> listOfStringTest();


    @SystemProperty
    @Key("system.Boolean.Def.True")
    @DefaultValue("true")
    boolean systemBooleanDefTrue();

    @SystemProperty
    @Key("system.Boolean.Def.False")
    @DefaultValue("false")
    boolean systemBooleanDefFalse();

    @SystemProperty
    @Key("system.Int.Def.207")
    @DefaultValue("207")
    int systemIntDef207();

    @SystemProperty
    @Key("system.List.Of.String.Test")
    @DefaultValue("string1,string2,string3,string4")
    List<String> systemListOfStringTest();
}
