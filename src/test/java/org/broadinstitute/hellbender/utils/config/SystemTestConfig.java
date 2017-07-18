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

    @DefaultValue("true")
    boolean systemBooleanDefTrue();

    @DefaultValue("false")
    boolean systemBooleanDefFalse();

    @DefaultValue("207")
    int systemIntDef207();

    @DefaultValue("string1,string2,string3,string4")
    List<String> systemListOfStringTest();



    @Key("system.Boolean.Def.True")
    @DefaultValue("true")
    boolean systemBooleanDefTrue2();

    @Key("system.Boolean.Def.False")
    @DefaultValue("false")
    boolean systemBooleanDefFalse2();

    @Key("system.Int.Def.207")
    @DefaultValue("207")
    int systemIntDef2072();

    @Key("system.List.Of.String.Test")
    @DefaultValue("string1,string2,string3,string4")
    List<String> systemListOfStringTest2();
}
