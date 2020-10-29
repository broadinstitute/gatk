{
<#--- This template is used only by tests, to generate test input data--->

<#--- Store positional args in a WDL arg called "positionalArgs"--->
<#assign positionalArgs="--positionalArgs"/>
  "${name}.dockerImage": "broadinstitute/gatk:${version}",
  "${name}.gatk": "java -cp ${buildDir}/build/libs/gatk.jar org.broadinstitute.hellbender.CommandLineArgumentValidatorMain",
<#if workflowProperties?? && workflowProperties?size != 0 && workflowProperties.memoryRequirements != "">
  "${name}.memoryRequirements": "${workflowProperties.memoryRequirements}",
<#else>
  "${name}.memoryRequirements": "String",
</#if>
<#if workflowProperties?? && workflowProperties?size != 0 && workflowProperties.diskRequirements != "">
  "${name}.diskRequirements": "${workflowProperties.diskRequirements}",
<#else>
  "${name}.diskRequirements": "String",
</#if>
<#if workflowProperties?? && workflowProperties?size != 0 && workflowProperties.cpuRequirements != "">
  "${name}.cpuRequirements": "${workflowProperties.cpuRequirements}",
<#else>
  "${name}.cpuRequirements": "String",
</#if>
<#if workflowProperties?? && workflowProperties?size != 0 && workflowProperties.preemptibleRequirements != "">
  "${name}.preemptibleRequirements": "${workflowProperties.preemptibleRequirements}",
<#else>
  "${name}.preemptibleRequirements": "String",
</#if>
<#if workflowProperties?? && workflowProperties?size != 0 && workflowProperties.bootdisksizegbRequirements != "">
  "${name}.bootdisksizegbRequirements": "${workflowProperties.bootdisksizegbRequirements}",
<#else>
  "${name}.bootdisksizegbRequirements": "String",
</#if>

<#assign remainingArgCount=arguments.required?size + arguments.optional?size + arguments.common?size/>
<@taskinput heading="Positional Arguments" argsToUse=arguments.positional remainingCount=remainingArgCount/>
<#assign remainingArgCount=arguments.optional?size + arguments.common?size/>
<@taskinput heading="Required Arguments" argsToUse=arguments.required remainingCount=remainingArgCount/>
<#assign remainingArgCount=arguments.common?size/>
<@taskinput heading="Optional Tool Arguments" argsToUse=arguments.optional remainingCount=remainingArgCount/>
<#assign remainingArgCount=0/>
<@taskinput heading="Optional Common Arguments" argsToUse=arguments.common remainingCount=0/>

}
<#macro taskinput heading argsToUse remainingCount>
  <#if argsToUse?size != 0>
    <#list argsToUse as arg>
      <#if heading?starts_with("Positional")>
          <#if requiredCompanions?? && requiredCompanions[positionalArgs]??>
              <#list requiredCompanions[positionalArgs] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: ${arg.testValue},
              </#list>
          </#if>
          <#if optionalCompanions?? && optionalCompanions[positionalArgs]??>
              <#list optionalCompanions[positionalArgs] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: ${arg.testValue},
              </#list>
          </#if>
<#noparse>  "</#noparse>${name}.${positionalArgs?substring(2)}<#noparse>"</#noparse>: <#rt/>
      <#else>
          <#if requiredCompanions?? && requiredCompanions[arg.name]??>
              <#list requiredCompanions[arg.name] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: ${arg.testValue},
              </#list>
          </#if>
          <#if optionalCompanions?? && optionalCompanions[arg.name]??>
              <#list optionalCompanions[arg.name] as companion>
<#noparse>  "</#noparse>${name}.${companion.name?substring(2)}<#noparse>"</#noparse>: ${arg.testValue},
              </#list>
          </#if>
<#noparse>  "</#noparse>${name}.${arg.name?substring(2)}<#noparse>"</#noparse>: <#rt/>
      </#if>
${arg.testValue}<#if !arg?is_last || remainingCount != 0>,</#if>
      <#if arg?is_last && remainingCount != 0>

      </#if>
    </#list>
</#if>
</#macro>
