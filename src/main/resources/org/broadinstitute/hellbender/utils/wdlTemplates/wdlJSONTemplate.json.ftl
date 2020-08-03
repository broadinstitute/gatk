{
<#--- Store positional args in a WDL arg called "positionalArgs"--->
<#assign positionalArgs="positionalArgs"/>
  "${name}.dockerImage": "broadinstitute/gatk:${version}",
  "${name}.gatk": "gatk",
<#if runtimeProperties?? && runtimeProperties?size != 0 && runtimeProperties.memoryRequirements != "">
  "${name}.memoryRequirements": "${runtimeProperties.memoryRequirements}",
<#else>
  "${name}.memoryRequirements": "String",
</#if>
<#if runtimeProperties?? && runtimeProperties?size != 0 && runtimeProperties.diskRequirements != "">
  "${name}.diskRequirements": "${runtimeProperties.diskRequirements}",
<#else>
  "${name}.diskRequirements": "String",
</#if>
<#if runtimeProperties?? && runtimeProperties?size != 0 && runtimeProperties.cpuRequirements != "">
  "${name}.cpuRequirements": "${runtimeProperties.cpuRequirements}",
<#else>
  "${name}.cpuRequirements": "String",
</#if>
<#if runtimeProperties?? && runtimeProperties?size != 0 && runtimeProperties.preemptibleRequirements != "">
  "${name}.preemptibleRequirements": "${runtimeProperties.preemptibleRequirements}",
<#else>
  "${name}.preemptibleRequirements": "String",
</#if>
<#if runtimeProperties?? && runtimeProperties?size != 0 && runtimeProperties.bootdisksizegbRequirements != "">
  "${name}.bootdisksizegbRequirements": "${runtimeProperties.bootdisksizegbRequirements}",
<#else>
  "${name}.bootdisksizegbRequirements": "String",
</#if>

<#assign remainingArgCount=arguments.required?size/>
<@taskinput heading="Positional Arguments" argsToUse=arguments.positional remainingCount=remainingArgCount/>
<#assign remainingArgCount=0/>
<@taskinput heading="Required Arguments" argsToUse=arguments.required remainingCount=remainingArgCount/>

}
<#macro taskinput heading argsToUse remainingCount>
  <#if argsToUse?size != 0>
    <#list argsToUse as arg>
      <#if heading?starts_with("Positional")>
          <#if requiredCompanions?? && requiredCompanions[positionalArgs]??>
              <#list requiredCompanions[positionalArgs] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: "${arg.wdlinputtype}",
              </#list>
          </#if>
          <#if optionalCompanions?? && optionalCompanions[positionalArgs]??>
              <#list optionalCompanions[positionalArgs] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: "${arg.wdlinputtype}",
              </#list>
          </#if>
<#noparse>  "</#noparse>${name}.${positionalArgs}<#noparse>"</#noparse>: <#rt/>
      <#else>
          <#if requiredCompanions?? && requiredCompanions[arg.name]??>
              <#list requiredCompanions[arg.name] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: "${arg.wdlinputtype}",
              </#list>
          </#if>
          <#if optionalCompanions?? && optionalCompanions[arg.name]??>
              <#list optionalCompanions[arg.name] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: "${arg.wdlinputtype}",
              </#list>
          </#if>
<#noparse>  "</#noparse>${name}.${arg.name?substring(2)}<#noparse>"</#noparse>: <#rt/>
      </#if>
<#noparse>  "</#noparse>${arg.wdlinputtype}<#noparse>"</#noparse><#if !arg?is_last || remainingCount != 0>,
    </#if>
      <#if arg?is_last && remainingCount != 0>

      </#if>
    </#list>
  </#if>
</#macro>
