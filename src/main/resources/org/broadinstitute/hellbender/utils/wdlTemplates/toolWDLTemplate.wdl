version 1.0

<#--- Store positional args in a WDL arg called "positionalArgs"--->
<#assign positionalArgs="positionalArgs"/>
<#if beta?? && beta == true>
# Run ${name} (**BETA**) (GATK Version ${version})
<#elseif experimental?? && experimental == true>
# Run ${name} **EXPERIMENTAL** ${name} (GATK Version ${version})
<#else>
# Run ${name} (GATK Version ${version})
</#if>
#
# ${summary}
#
# General Workflow (non-tool) Arguments
#    ${"gatk"?right_pad(50)}Location of gatk used to run this workflow
#
<#if arguments.positional?size != 0>
<@addArgumentDescriptions heading="Positional Tool Arguments" argsToUse=arguments.positional/>
#
</#if>
<#if arguments.required?size != 0>
<@addArgumentDescriptions heading="Required Tool Arguments" argsToUse=arguments.required/>
#
</#if>
<#if arguments.optional?size != 0>
<@addArgumentDescriptions heading="Optional Tool Arguments" argsToUse=arguments.optional/>
#
</#if>
<#if arguments.common?size != 0>
<@addArgumentDescriptions heading="Optional Common Arguments" argsToUse=arguments.common/>
</#if>

workflow ${name} {

  input {
    #GATK location
    String gatk
    <@defineWorkflowInputs heading="Positional Arguments" argsToUse=arguments.positional/>
    <@defineWorkflowInputs heading="Required Arguments" argsToUse=arguments.required/>
    <@defineWorkflowInputs heading="Optional Tool Arguments" argsToUse=arguments.optional/>
    <@defineWorkflowInputs heading="Optional Common Arguments" argsToUse=arguments.common/>

  }

  call ${name}Task {

    input:

        #GATK location
        ${"gatk"?right_pad(50)} = gatk,
        <@callTaskInputs heading="Positional Arguments" argsToUse=arguments.positional/>
        <@callTaskInputs heading="Required Arguments" argsToUse=arguments.required/>
        <@callTaskInputs heading="Optional Tool Arguments" argsToUse=arguments.optional/>
        <@callTaskInputs heading="Optional Common Arguments" argsToUse=arguments.common/>

  }

  output {
    <@defineWorkflowOutputs heading="Workflow Outputs" outputs=RuntimeOutputs/>
  }
}

task ${name}Task {

  input {
    String gatk
    <@defineTaskInputs heading="Positional Arguments" argsToUse=arguments.positional/>
    <@defineTaskInputs heading="Required Arguments" argsToUse=arguments.required/>
    <@defineTaskInputs heading="Optional Tool Arguments" argsToUse=arguments.optional/>
    <@defineTaskInputs heading="Optional Common Arguments" argsToUse=arguments.common/>

  }

  command <<<
    ~{gatk} ${name} \
        <@callTaskCommand heading="Positional Arguments" argsToUse=arguments.positional/>
        <@callTaskCommand heading="Required Arguments" argsToUse=arguments.required/>
        <@callTaskCommand heading="Optional Tool Arguments" argsToUse=arguments.optional/>
        <@callTaskCommand heading="Optional Common Arguments" argsToUse=arguments.common/>
  >>>

  <#if runtimeProperties?? && runtimeProperties?size != 0>
  runtime {
      <#if runtimeProperties.memory != "">
    memory: "${runtimeProperties.memory}"
      </#if>
  }
  </#if>

  output {
    <@defineTaskOutputs heading="Task Outputs" outputs=RuntimeOutputs/>
  }
 }

<#--------------------------------------->
<#-- Macros -->

<#macro addArgumentDescriptions heading argsToUse>
    <#if argsToUse?size != 0>
# ${heading}
        <#list argsToUse as arg>
#   ${arg.name?substring(2)?right_pad(50)} ${arg.summary?right_pad(60)[0..*80]}
        </#list>
    </#if>
</#macro>

<#macro defineWorkflowInputs heading argsToUse>
    <#if argsToUse?size != 0>

    # ${heading}
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
    ${arg.type} ${positionalArgs}
            <#else>
    ${arg.type}<#if !heading?starts_with("Required")>?</#if> ${arg.name?substring(2)}
            </#if>
        </#list>
    </#if>
</#macro>

<#macro callTaskInputs heading argsToUse>
    <#if argsToUse?size != 0>

        # ${heading}
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
        ${positionalArgs?right_pad(50)} = ${positionalArgs},
            <#else>
        ${arg.name?substring(2)?right_pad(50)} = ${arg.name?substring(2)},
            </#if>
        </#list>
    </#if>
</#macro>

<#macro defineTaskInputs heading argsToUse>
    <#if argsToUse?size != 0>
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
    ${arg.type} positionalArgs
            <#else>
    ${arg.type}<#if !heading?starts_with("Required")>?</#if> ${arg.name?substring(2)}
            </#if>
        </#list>
    </#if>
</#macro>

<#macro defineWorkflowOutputs heading outputs>
    # ${heading?right_pad(50)}
    <#if outputs?size == 0>
    File ${name}results = ${name}Task.${name}Task_results
    <#else>
        <#list outputs as outputFile>
    File ${name}${outputFile} = ${name}Task.${name}Task_${outputFile}
        </#list>
    </#if>
</#macro>

<#macro defineTaskOutputs heading outputs>
    # ${heading?right_pad(50)}
    <#if outputs?size == 0>
    File ${name}Task_results = stdout()
    <#else>
        <#list outputs as outputFile>
    File ${name}Task_${outputFile} = <#noparse>"${</#noparse>${outputFile}<#noparse>}"</#noparse>
        </#list>
    </#if>
</#macro>

<#macro callTaskCommand heading argsToUse>
    <#if argsToUse?size != 0>
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
        <#noparse>~{sep=' ' </#noparse>${positionalArgs}<#noparse>}</#noparse> \
            <#elseif heading?starts_with("Required")>
        ${arg.name} <#noparse>~{sep=' </#noparse>${arg.name} <#noparse>' </#noparse>${arg.name?substring(2)}<#noparse>}</#noparse> \
            <#else>
        <#noparse>~{true='</#noparse>${arg.name} <#noparse>' false='' defined(</#noparse>${arg.name?substring(2)}<#noparse>)}~{sep='</#noparse> ${arg.name} <#noparse>'</#noparse> ${arg.name?substring(2)}<#noparse>}</#noparse> \
            </#if>
        </#list>
    </#if>
</#macro>
