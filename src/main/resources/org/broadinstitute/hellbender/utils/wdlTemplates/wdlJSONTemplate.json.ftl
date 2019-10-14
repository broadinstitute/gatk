{
<#--- Store positional args in a WDL arg called "positionalArgs"--->
<#assign positionalArgs="positionalArgs"/>
  "${name}.gatk": "String",
<#assign remainingArgCount=arguments.required?size + arguments.optional?size + arguments.common?size/>
<@taskinput heading="Positional Arguments" argsToUse=arguments.positional remainingCount=remainingArgCount/>
<#assign remainingArgCount=arguments.optional?size + arguments.common?size/>
<@taskinput heading="Required Arguments" argsToUse=arguments.required remainingCount=remainingArgCount/>
<#assign remainingArgCount=arguments.common?size/>
<@taskinput heading="Optional Tool Arguments" argsToUse=arguments.optional remainingCount=remainingArgCount/>
<#assign remainingArgCount=arguments.required?size + arguments.optional?size + arguments.common?size/>
<@taskinput heading="Optional Common Arguments" argsToUse=arguments.common remainingCount=0/>
}
<#macro taskinput heading argsToUse remainingCount>
  <#if argsToUse?size != 0>
    <#list argsToUse as arg>
      <#if companionResources?? && companionResources[arg.name]??>
            <#list companionResources[arg.name] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: null,
            </#list>
        </#if>
    <#if heading?starts_with("Positional")>
<#noparse>  "</#noparse>${name}.${positionalArgs}<#noparse>"</#noparse>: <#rt/>
      <#else>
<#noparse>  "</#noparse>${name}.${arg.name?substring(2)}<#noparse>"</#noparse>: <#rt/>
      </#if>
      <#if heading?starts_with("Required") || heading?starts_with("Positional")>
<#noparse>"</#noparse>${arg.type}<#noparse>"</#noparse><#if !arg?is_last || remainingCount != 0>,</#if>
      <#else>
        <#if arg.defaultValue == "[]" || arg.defaultValue == "\"\"" || arg.defaultValue == "null">
null<#if !arg?is_last || remainingCount != 0>,</#if>
        <#else>
<#noparse>"</#noparse>${arg.defaultValue}<#noparse>"</#noparse><#if !arg?is_last || remainingCount != 0>,</#if>
        </#if>
      </#if>
      <#if arg?is_last && remainingCount != 0>

      </#if>
    </#list>
  </#if>
</#macro>
