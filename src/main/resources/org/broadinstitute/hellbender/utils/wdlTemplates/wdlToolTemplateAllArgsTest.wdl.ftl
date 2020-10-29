<#--- This template is used only by tests, and differs from the *AllArgs WDL--->
<#--- in that is has no runtime outputs or runtime docker declared --->
version 1.0

<#--- Store positional args in a WDL arg called "positionalArgs"--->
<#assign positionalArgs="--positionalArgs"/>
<#if beta?? && beta == true>
# Run ${name} (**BETA**) (WDL auto generated from: GATK Version ${version})
<#elseif experimental?? && experimental == true>
# Run ${name} **EXPERIMENTAL** ${name} (WDL auto generated from: GATK Version ${version})
<#else>
# Run ${name} (WDL auto generated from: GATK Version ${version})
</#if>
#
# ${summary}
#
#  General Workflow (non-tool) Arguments
#    ${"dockerImage"?right_pad(50)} Docker image for this workflow
#    ${"gatk"?right_pad(50)} Location of gatk to run for this workflow
#    ${"memoryRequirements"?right_pad(50)} Runtime memory requirements for this workflow
#    ${"diskRequirements"?right_pad(50)} Runtime disk requirements for this workflow
#    ${"cpuRequirements"?right_pad(50)} Runtime CPU count for this workflow
#    ${"preemptibleRequirements"?right_pad(50)} Runtime preemptible count for this workflow
#    ${"bootdisksizegbRequirements"?right_pad(50)} Runtime boot disk size for this workflow
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
    #Docker to use
    String dockerImage
    #App location
    String gatk
    #Memory to use
    String memoryRequirements
    #Disk requirements for this workflow
    String diskRequirements
    #CPU requirements for this workflow
    String cpuRequirements
    #Preemptible requirements for this workflow
    String preemptibleRequirements
    #Boot disk size requirements for this workflow
    String bootdisksizegbRequirements
    <@defineWorkflowInputs heading="Positional Arguments" argsToUse=arguments.positional/>
    <@defineWorkflowInputs heading="Required Arguments" argsToUse=arguments.required/>
    <@defineWorkflowInputs heading="Optional Tool Arguments" argsToUse=arguments.optional/>
    <@defineWorkflowInputs heading="Optional Common Arguments" argsToUse=arguments.common/>

  }

  call ${name} {

    input:

        #Docker
        ${"dockerImage"?right_pad(50)} = dockerImage,
        #App location
        ${"gatk"?right_pad(50)} = gatk,
        #Memory to use
        ${"memoryRequirements"?right_pad(50)} = memoryRequirements,
        #Disk requirements for this workflow
        ${"diskRequirements"?right_pad(50)} = diskRequirements,
        #CPU requirements for this workflow
        ${"cpuRequirements"?right_pad(50)} = cpuRequirements,
        #Preemptible requirements for this workflow
        ${"preemptibleRequirements"?right_pad(50)} = preemptibleRequirements,
        #Boot disk size requirements for this workflow
        ${"bootdisksizegbRequirements"?right_pad(50)} = bootdisksizegbRequirements,

        <@callTaskInputs heading="Positional Arguments" argsToUse=arguments.positional/>
        <@callTaskInputs heading="Required Arguments" argsToUse=arguments.required/>
        <@callTaskInputs heading="Optional Tool Arguments" argsToUse=arguments.optional/>
        <@callTaskInputs heading="Optional Common Arguments" argsToUse=arguments.common/>

  }

  parameter_meta {
    dockerImage: { description: "Docker image for this task" }
    gatk: { description: "Location of gatk to run for this task" }
    memoryRequirements: { description: "Runtime memory requirements for this task" }
    diskRequirements: { description: "Runtime disk requirements for this task" }
    cpuRequirements: { description: "Runtime CPU count for this task" }
    preemptibleRequirements: { description: "Runtime preemptible count for this task" }
    bootdisksizegbRequirements: { description: "Runtime boot disk size for this task" }
    <@defineParamMeta heading="Positional Arguments" argsToUse=arguments.positional/>
    <@defineParamMeta heading="Required Arguments" argsToUse=arguments.required/>
    <@defineParamMeta heading="Optional Tool Arguments" argsToUse=arguments.optional/>
    <@defineParamMeta heading="Optional Common Arguments" argsToUse=arguments.common/>
  }
}

task ${name} {

  input {
    String dockerImage
    String gatk
    String memoryRequirements
    String diskRequirements
    String cpuRequirements
    String preemptibleRequirements
    String bootdisksizegbRequirements
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

  <#if workflowProperties?? && workflowProperties?size != 0>
  runtime {
      #docker: dockerImage
      memory: memoryRequirements
      disks: diskRequirements
      cpu: cpuRequirements
      preemptible: preemptibleRequirements
      bootDiskSizeGb: bootdisksizegbRequirements
  }
  </#if>

  parameter_meta {
    dockerImage: { description: "Docker image for this task" }
    gatk: { description: "Location of gatk to run for this task" }
    memoryRequirements: { description: "Runtime memory requirements for this task" }
    diskRequirements: { description: "Runtime disk requirements for this task" }
    cpuRequirements: { description: "Runtime CPU count for this task" }
    preemptibleRequirements: { description: "Runtime preemptible count for this task" }
    bootdisksizegbRequirements: { description: "Runtime boot disk size for this task" }
    <@defineParamMeta heading="Positional Arguments" argsToUse=arguments.positional/>
    <@defineParamMeta heading="Required Arguments" argsToUse=arguments.required/>
    <@defineParamMeta heading="Optional Tool Arguments" argsToUse=arguments.optional/>
    <@defineParamMeta heading="Optional Common Arguments" argsToUse=arguments.common/>
  }
}

<#--------------------------------------->
<#-- Macros -->

<#macro addArgumentDescriptions heading argsToUse>
    <#if argsToUse?size != 0>
#  ${heading}
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
#    ${positionalArgs?substring(2)?right_pad(50)} ${arg.summary?right_pad(60)[0..*80]}
                <#if requiredCompanions?? && requiredCompanions[positionalArgs]??>
                    <#list requiredCompanions[positionalArgs] as companion>
#    ${companion.name?substring(2)?right_pad(50)} ${companion.summary?right_pad(60)[0..*80]}
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[positionalArgs]??>
                    <#list optionalCompanions[positionalArgs] as companion>
#    ${companion.name?substring(2)?right_pad(50)} Optional ${companion.summary?right_pad(60)[0..*80]}
                    </#list>
                </#if>
            <#else>
#    ${arg.name?substring(2)?right_pad(50)} ${arg.summary?right_pad(60)[0..*80]}
                <#if requiredCompanions?? && requiredCompanions[arg.name]??>
                    <#list requiredCompanions[arg.name] as companion>
#    ${companion.name?substring(2)?right_pad(50)} ${companion.summary?right_pad(60)[0..*80]}
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[arg.name]??>
                    <#list optionalCompanions[arg.name] as companion>
#    ${companion.name?substring(2)?right_pad(50)} Optional ${companion.summary?right_pad(60)[0..*80]}
                    </#list>
                </#if>
            </#if>
        </#list>
    </#if>
</#macro>


<#macro defineParamMeta heading argsToUse>
    <#if argsToUse?size != 0>

    # ${heading}
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
                <#if localizationOptional?seq_contains(positionalArgs)>
    ${positionalArgs?substring(2)}<#noparse>: {
      description: "</#noparse>${arg.summary[0..*80]}<#noparse>",</#noparse>
      localization_optional : true <#noparse>
    }</#noparse>
                <#else>
    ${positionalArgs?substring(2)}<#noparse>: { description: "</#noparse>${arg.summary[0..*80]}<#noparse>" }</#noparse>
                </#if>
                <#if requiredCompanions?? && requiredCompanions[positionalArgs]??>
                    <#list requiredCompanions[positionalArgs] as companion>
                        <#if localizationOptional?seq_contains(positionalArgs)>
    ${companion.name?substring(2)}<#noparse>: {
      description: "</#noparse>${companion.summary[0..*80]}<#noparse>",</#noparse>
      localization_optional : true <#noparse>
    }</#noparse>
                        <#else>
    ${companion.name?substring(2)}<#noparse>: { description: "</#noparse>${companion.summary[0..*80]}<#noparse>" }</#noparse>
                        </#if>
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[positionalArgs]??>
                    <#list optionalCompanions[positionalArgs] as companion>
                        <#if localizationOptional?seq_contains(positionalArgs)>
    ${companion.name?substring(2)}<#noparse>: {
      description: "</#noparse>${arg.summary[0..*80]}<#noparse>",</#noparse>
      localization_optional : true <#noparse>
    }</#noparse>
                        <#else>
    ${companion.name?substring(2)}<#noparse>: { description: "</#noparse>${companion.summary(60)[0..*80]}<#noparse>" }</#noparse>
                        </#if>
                    </#list>
                </#if>
            <#else>
                <#if localizationOptional?seq_contains(arg.name)>
    ${arg.name?substring(2)}<#noparse>: {
      description: "</#noparse>${arg.summary[0..*80]}<#noparse>",</#noparse>
      localization_optional : true <#noparse>
    }</#noparse>
                <#else>
    ${arg.name?substring(2)}<#noparse>: { description: "</#noparse>${arg.summary[0..*80]}<#noparse>" }</#noparse>
                </#if>
                <#if requiredCompanions?? && requiredCompanions[arg.name]??>
                    <#list requiredCompanions[arg.name] as companion>
                        <#if localizationOptional?seq_contains(companion.name)>
    ${companion.name?substring(2)}<#noparse>: {
      description: "</#noparse>${companion.summary[0..*80]}<#noparse>",</#noparse>
      localization_optional : true <#noparse>
    }</#noparse>
                        <#else>
    ${companion.name?substring(2)}<#noparse>: { description: "</#noparse>${companion.summary[0..*80]}<#noparse>" }</#noparse>
                        </#if>
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[arg.name]??>
                    <#list optionalCompanions[arg.name] as companion>
                        <#if localizationOptional?seq_contains(companion.name)>
    ${companion.name?substring(2)}<#noparse>: {
      description: "</#noparse>${companion.summary[0..*80]}<#noparse>",</#noparse>
      localization_optional : true <#noparse>
    }</#noparse>
                        <#else>
    ${companion.name?substring(2)}<#noparse>: { description: "</#noparse>${companion.summary[0..*80]}<#noparse>" }</#noparse>
                        </#if>
                    </#list>
                </#if>
            </#if>
        </#list>
    </#if>
</#macro>

<#macro defineWorkflowInputs heading argsToUse>
    <#if argsToUse?size != 0>

    # ${heading}
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
    ${arg.wdlinputtype} ${positionalArgs?substring(2)}
                <#if requiredCompanions?? && requiredCompanions[positionalArgs]??>
                    <#list requiredCompanions[positionalArgs] as companion>
    ${arg.wdlinputtype} ${companion.name?substring(2)}
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[positionalArgs]??>
                    <#list optionalCompanions[positionalArgs] as companion>
    ${arg.wdlinputtype}? ${companion.name?substring(2)}
                    </#list>
                </#if>
            <#else>
    ${arg.wdlinputtype}<#if !heading?starts_with("Required")>?</#if> ${arg.name?substring(2)}
                <#if requiredCompanions?? && requiredCompanions[arg.name]??>
                    <#list requiredCompanions[arg.name] as companion>
    ${arg.wdlinputtype} ${companion.name?substring(2)}
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[arg.name]??>
                    <#list optionalCompanions[arg.name] as companion>
    ${arg.wdlinputtype}? ${companion.name?substring(2)}
                    </#list>
                </#if>
            </#if>
        </#list>
    </#if>
</#macro>

<#macro callTaskInputs heading argsToUse>
    <#if argsToUse?size != 0>

        # ${heading}
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
        ${positionalArgs?substring(2)?right_pad(50)} = ${positionalArgs?substring(2)},
                <#if requiredCompanions?? && requiredCompanions[positionalArgs]??>
                    <#list requiredCompanions[positionalArgs] as companion>
        ${companion.name?substring(2)?right_pad(50)} = ${companion.name?substring(2)},
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[positionalArgs]??>
                    <#list optionalCompanions[positionalArgs] as companion>
        ${companion.name?substring(2)?right_pad(50)} = ${companion.name?substring(2)},
                    </#list>
                </#if>
            <#else>
        ${arg.name?substring(2)?right_pad(50)} = ${arg.name?substring(2)},
                <#if requiredCompanions?? && requiredCompanions[arg.name]??>
                    <#list requiredCompanions[arg.name] as companion>
        ${companion.name?substring(2)?right_pad(50)} = ${companion.name?substring(2)},
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[arg.name]??>
                    <#list optionalCompanions[arg.name] as companion>
        ${companion.name?substring(2)?right_pad(50)} = ${companion.name?substring(2)},
                    </#list>
                </#if>
            </#if>
        </#list>
    </#if>
</#macro>

<#macro defineTaskInputs heading argsToUse>
    <#if argsToUse?size != 0>
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
    ${arg.wdlinputtype} ${positionalArgs?substring(2)}
                <#if requiredCompanions?? && requiredCompanions[arg.name]??>
                    <#list requiredCompanions[positionalArgs] as companion>
    ${arg.wdlinputtype} Positional_${companion.name?substring(2)}
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[arg.name]??>
                    <#list optionalCompanions[positionalArgs] as companion>
    ${arg.wdlinputtype}? Positional_${companion.name?substring(2)}
                    </#list>
                </#if>
            <#else>
    ${arg.wdlinputtype}<#if !heading?starts_with("Required")>?</#if> ${arg.name?substring(2)}
                <#if requiredCompanions?? && requiredCompanions[arg.name]??>
                    <#list requiredCompanions[arg.name] as companion>
    ${arg.wdlinputtype} ${companion.name?substring(2)}
                    </#list>
                </#if>
                <#if optionalCompanions?? && optionalCompanions[arg.name]??>
                    <#list optionalCompanions[arg.name] as companion>
    ${arg.wdlinputtype}? ${companion.name?substring(2)}
                    </#list>
                </#if>
            </#if>
        </#list>
    </#if>
</#macro>

<#macro defineWorkflowOutputs heading outputs>
    # ${heading?right_pad(50)}
    <#if outputs?size == 0>
    File ${name}results = ${name}.${name}_results
    <#else>
        <#list outputs?keys as key>
            <#if requiredOutputs[key]??>
    ${outputs[key]} ${name}${key?substring(2)} = ${name}.${name}_${key?substring(2)}
            <#else>
    ${outputs[key]}? ${name}${key?substring(2)} = ${name}.${name}_${key?substring(2)}
            </#if>
        </#list>
    </#if>
</#macro>

<#macro defineTaskOutputs heading outputs>
    # ${heading?right_pad(50)}
    <#if outputs?size == 0>
    File ${name}_results = stdout()
    <#else>
        <#list outputs?keys as key>
            <#if requiredOutputs[key]??>
    ${outputs[key]} ${name}_${key?substring(2)} = ${key?substring(2)}
            <#else>
    ${outputs[key]}? ${name}_${key?substring(2)} = ${key?substring(2)}
            </#if>
        </#list>
    </#if>
</#macro>

<#macro callTaskCommand heading argsToUse>
    <#if argsToUse?size != 0>
        <#list argsToUse as arg>
            <#if heading?starts_with("Positional")>
    <#noparse>~{sep=' ' </#noparse>${positionalArgs?substring(2)}<#noparse>}</#noparse> \
            <#elseif heading?starts_with("Required")>
    ${arg.actualArgName} <#noparse>~{sep=' </#noparse>${arg.actualArgName} <#noparse>' </#noparse>${arg.name?substring(2)}<#noparse>}</#noparse> \
            <#else>
    <#noparse>~{true='</#noparse>${arg.actualArgName} <#noparse>' false='' defined(</#noparse>${arg.name?substring(2)}<#noparse>)}~{sep='</#noparse> ${arg.actualArgName} <#noparse>'</#noparse> ${arg.name?substring(2)}<#noparse>}</#noparse> \
            </#if>
        </#list>
    </#if>
</#macro>
