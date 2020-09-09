version 1.0

# Run TestWDLTool (WDL auto generated from GATK Version 1.2.3.4)
#
# WDL Test Tool to test WDL Generation
#
#  General Workflow (non-tool) Arguments
#    dockerImage                                        Docker image for this workflow
#    gatk                                               Location of gatk to run for this workflow
#    memoryRequirements                                 Runtime memory requirements for this workflow
#    diskRequirements                                   Runtime disk requirements for this workflow
#    cpuRequirements                                    Runtime CPU count for this workflow
#    preemptibleRequirements                            Runtime preemptible count for this workflow
#    bootdisksizegbRequirements                         Runtime boot disk size for this workflow
#
#  Positional Tool Arguments
#    positionalArgs                                     Positional args doc                                         
#    posDictionary                                      Companion resource for positionalArgs                       
#    posIndex                                           Companion resource for positionalArgs                       
#
#  Required Tool Arguments
#    requiredListFileInputMixedCompanions               requiredListFileInputMixedCompanions doc                    
#    requiredListFileInputMixedCompanionsRequired       Companion resource for requiredListFileInputMixedCompanions 
#    requiredListFileInputMixedCompanionsOptional       Optional Companion resource for requiredListFileInputMixedCompanions 
#    requiredListFileInputNoCompanions                  requiredListFileInputNoCompanions doc                       
#    requiredListFileInputOptionalCompanions            requiredListFileInputOptionalCompanions doc                 
#    requiredListFileInputOptionalCompanionsDictionary  Optional Companion resource for requiredListFileInputOptionalCompanions
#    requiredListFileInputOptionalCompanionsIndex       Optional Companion resource for requiredListFileInputOptionalCompanions
#    requiredListFileInputRequiredCompanions            requiredListFileInputRequiredCompanions doc                 
#    requiredListFileInputRequiredCompanionsDictionary  Companion resource for requiredListFileInputRequiredCompanions
#    requiredListFileInputRequiredCompanionsIndex       Companion resource for requiredListFileInputRequiredCompanions
#    requiredListFileOutputMixedCompanions              requiredListFileOutputMixedCompanions doc                   
#    requiredListFileOutputMixedCompanionsRequired      Companion resource for requiredListFileOutputMixedCompanions
#    requiredListFileOutputMixedCompanionsOptional      Optional Companion resource for requiredListFileOutputMixedCompanions
#    requiredListFileOutputNoCompanions                 requiredListFileOutputNoCompanions doc                      
#    requiredListFileOutputOptionalCompanions           requiredListFileOutputOptionalCompanions doc                
#    requiredListFileOutputOptionalCompanionsDictionary Optional Companion resource for requiredListFileOutputOptionalCompanions
#    requiredListFileOutputOptionalCompanionsIndex      Optional Companion resource for requiredListFileOutputOptionalCompanions
#    requiredListFileOutputRequiredCompanions           requiredListFileOutputRequiredCompanions doc                
#    requiredListFileOutputRequiredCompanionsDictionary Companion resource for requiredListFileOutputRequiredCompanions
#    requiredListFileOutputRequiredCompanionsIndex      Companion resource for requiredListFileOutputRequiredCompanions
#    requiredScalarFileInputNoCompanions                requiredScalarFileInputNoCompanions doc                     
#    requiredScalarFileInputOptionalCompanions          requiredScalarFileInputOptionalCompanions doc               
#    requiredScalarFileInputOptionalCompanionsDictionary Optional Companion resource for requiredScalarFileInputOptionalCompanions
#    requiredScalarFileInputOptionalCompanionsIndex     Optional Companion resource for requiredScalarFileInputOptionalCompanions
#    requiredScalarFileInputRequiredCompanions          requiredScalarFileInputRequiredCompanions doc               
#    requiredScalarFileInputRequiredCompanionsDictionary Companion resource for requiredScalarFileInputRequiredCompanions
#    requiredScalarFileInputRequiredCompanionsIndex     Companion resource for requiredScalarFileInputRequiredCompanions
#    requiredScalarFileOutputNoCompanions               requiredScalarFileOutputNoCompanions doc                    
#    requiredScalarFileOutputOptionalCompanions         requiredScalarFileOutputOptionalCompanions doc              
#    requiredScalarFileOutputOptionalCompanionsDictionary Optional Companion resource for requiredScalarFileOutputOptionalCompanions
#    requiredScalarFileOutputOptionalCompanionsIndex    Optional Companion resource for requiredScalarFileOutputOptionalCompanions
#    requiredScalarFileOutputRequiredCompanions         requiredScalarFileOutputRequiredCompanions doc              
#    requiredScalarFileOutputRequiredCompanionsDictionary Companion resource for requiredScalarFileOutputRequiredCompanions
#    requiredScalarFileOutputRequiredCompanionsIndex    Companion resource for requiredScalarFileOutputRequiredCompanions
#

workflow TestWDLTool {

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

    # Positional Arguments
    Array[File] positionalArgs
    Array[File] posDictionary
    Array[File] posIndex

    # Required Arguments
    Array[File] requiredListFileInputMixedCompanions
    Array[File] requiredListFileInputMixedCompanionsRequired
    Array[File]? requiredListFileInputMixedCompanionsOptional
    Array[File] requiredListFileInputNoCompanions
    Array[File] requiredListFileInputOptionalCompanions
    Array[File]? requiredListFileInputOptionalCompanionsDictionary
    Array[File]? requiredListFileInputOptionalCompanionsIndex
    Array[File] requiredListFileInputRequiredCompanions
    Array[File] requiredListFileInputRequiredCompanionsDictionary
    Array[File] requiredListFileInputRequiredCompanionsIndex
    Array[String] requiredListFileOutputMixedCompanions
    Array[String] requiredListFileOutputMixedCompanionsRequired
    Array[String]? requiredListFileOutputMixedCompanionsOptional
    Array[String] requiredListFileOutputNoCompanions
    Array[String] requiredListFileOutputOptionalCompanions
    Array[String]? requiredListFileOutputOptionalCompanionsDictionary
    Array[String]? requiredListFileOutputOptionalCompanionsIndex
    Array[String] requiredListFileOutputRequiredCompanions
    Array[String] requiredListFileOutputRequiredCompanionsDictionary
    Array[String] requiredListFileOutputRequiredCompanionsIndex
    File requiredScalarFileInputNoCompanions
    File requiredScalarFileInputOptionalCompanions
    File? requiredScalarFileInputOptionalCompanionsDictionary
    File? requiredScalarFileInputOptionalCompanionsIndex
    File requiredScalarFileInputRequiredCompanions
    File requiredScalarFileInputRequiredCompanionsDictionary
    File requiredScalarFileInputRequiredCompanionsIndex
    String requiredScalarFileOutputNoCompanions
    String requiredScalarFileOutputOptionalCompanions
    String? requiredScalarFileOutputOptionalCompanionsDictionary
    String? requiredScalarFileOutputOptionalCompanionsIndex
    String requiredScalarFileOutputRequiredCompanions
    String requiredScalarFileOutputRequiredCompanionsDictionary
    String requiredScalarFileOutputRequiredCompanionsIndex

  }

  call TestWDLTool {

    input:

        #Docker
        dockerImage                                        = dockerImage,
        #App location
        gatk                                               = gatk,
        #Memory to use
        memoryRequirements                                 = memoryRequirements,
        #Disk requirements for this workflow
        diskRequirements                                   = diskRequirements,
        #CPU requirements for this workflow
        cpuRequirements                                    = cpuRequirements,
        #Preemptible requirements for this workflow
        preemptibleRequirements                            = preemptibleRequirements,
        #Boot disk size requirements for this workflow
        bootdisksizegbRequirements                         = bootdisksizegbRequirements,


        # Positional Arguments
        positionalArgs                                     = positionalArgs,
        posDictionary                                      = posDictionary,
        posIndex                                           = posIndex,

        # Required Arguments
        requiredListFileInputMixedCompanions               = requiredListFileInputMixedCompanions,
        requiredListFileInputMixedCompanionsRequired       = requiredListFileInputMixedCompanionsRequired,
        requiredListFileInputMixedCompanionsOptional       = requiredListFileInputMixedCompanionsOptional,
        requiredListFileInputNoCompanions                  = requiredListFileInputNoCompanions,
        requiredListFileInputOptionalCompanions            = requiredListFileInputOptionalCompanions,
        requiredListFileInputOptionalCompanionsDictionary  = requiredListFileInputOptionalCompanionsDictionary,
        requiredListFileInputOptionalCompanionsIndex       = requiredListFileInputOptionalCompanionsIndex,
        requiredListFileInputRequiredCompanions            = requiredListFileInputRequiredCompanions,
        requiredListFileInputRequiredCompanionsDictionary  = requiredListFileInputRequiredCompanionsDictionary,
        requiredListFileInputRequiredCompanionsIndex       = requiredListFileInputRequiredCompanionsIndex,
        requiredListFileOutputMixedCompanions              = requiredListFileOutputMixedCompanions,
        requiredListFileOutputMixedCompanionsRequired      = requiredListFileOutputMixedCompanionsRequired,
        requiredListFileOutputMixedCompanionsOptional      = requiredListFileOutputMixedCompanionsOptional,
        requiredListFileOutputNoCompanions                 = requiredListFileOutputNoCompanions,
        requiredListFileOutputOptionalCompanions           = requiredListFileOutputOptionalCompanions,
        requiredListFileOutputOptionalCompanionsDictionary = requiredListFileOutputOptionalCompanionsDictionary,
        requiredListFileOutputOptionalCompanionsIndex      = requiredListFileOutputOptionalCompanionsIndex,
        requiredListFileOutputRequiredCompanions           = requiredListFileOutputRequiredCompanions,
        requiredListFileOutputRequiredCompanionsDictionary = requiredListFileOutputRequiredCompanionsDictionary,
        requiredListFileOutputRequiredCompanionsIndex      = requiredListFileOutputRequiredCompanionsIndex,
        requiredScalarFileInputNoCompanions                = requiredScalarFileInputNoCompanions,
        requiredScalarFileInputOptionalCompanions          = requiredScalarFileInputOptionalCompanions,
        requiredScalarFileInputOptionalCompanionsDictionary = requiredScalarFileInputOptionalCompanionsDictionary,
        requiredScalarFileInputOptionalCompanionsIndex     = requiredScalarFileInputOptionalCompanionsIndex,
        requiredScalarFileInputRequiredCompanions          = requiredScalarFileInputRequiredCompanions,
        requiredScalarFileInputRequiredCompanionsDictionary = requiredScalarFileInputRequiredCompanionsDictionary,
        requiredScalarFileInputRequiredCompanionsIndex     = requiredScalarFileInputRequiredCompanionsIndex,
        requiredScalarFileOutputNoCompanions               = requiredScalarFileOutputNoCompanions,
        requiredScalarFileOutputOptionalCompanions         = requiredScalarFileOutputOptionalCompanions,
        requiredScalarFileOutputOptionalCompanionsDictionary = requiredScalarFileOutputOptionalCompanionsDictionary,
        requiredScalarFileOutputOptionalCompanionsIndex    = requiredScalarFileOutputOptionalCompanionsIndex,
        requiredScalarFileOutputRequiredCompanions         = requiredScalarFileOutputRequiredCompanions,
        requiredScalarFileOutputRequiredCompanionsDictionary = requiredScalarFileOutputRequiredCompanionsDictionary,
        requiredScalarFileOutputRequiredCompanionsIndex    = requiredScalarFileOutputRequiredCompanionsIndex,

  }

  output {
    # Workflow Outputs                                  
    File TestWDLToolrequiredScalarFileOutputNoCompanions = TestWDLTool.TestWDLTool_requiredScalarFileOutputNoCompanions
    File TestWDLToolrequiredScalarFileOutputRequiredCompanions = TestWDLTool.TestWDLTool_requiredScalarFileOutputRequiredCompanions
    File TestWDLToolrequiredScalarFileOutputRequiredCompanionsDictionary = TestWDLTool.TestWDLTool_requiredScalarFileOutputRequiredCompanionsDictionary
    File TestWDLToolrequiredScalarFileOutputRequiredCompanionsIndex = TestWDLTool.TestWDLTool_requiredScalarFileOutputRequiredCompanionsIndex
    File TestWDLToolrequiredScalarFileOutputRequiredCompanionsDictionary = TestWDLTool.TestWDLTool_requiredScalarFileOutputRequiredCompanionsDictionary
    File TestWDLToolrequiredScalarFileOutputRequiredCompanionsIndex = TestWDLTool.TestWDLTool_requiredScalarFileOutputRequiredCompanionsIndex
    File TestWDLToolrequiredScalarFileOutputOptionalCompanions = TestWDLTool.TestWDLTool_requiredScalarFileOutputOptionalCompanions
    File? TestWDLToolrequiredScalarFileOutputOptionalCompanionsDictionary = TestWDLTool.TestWDLTool_requiredScalarFileOutputOptionalCompanionsDictionary
    File? TestWDLToolrequiredScalarFileOutputOptionalCompanionsIndex = TestWDLTool.TestWDLTool_requiredScalarFileOutputOptionalCompanionsIndex
    Array[File] TestWDLToolrequiredListFileOutputNoCompanions = TestWDLTool.TestWDLTool_requiredListFileOutputNoCompanions
    Array[File] TestWDLToolrequiredListFileOutputRequiredCompanions = TestWDLTool.TestWDLTool_requiredListFileOutputRequiredCompanions
    Array[File] TestWDLToolrequiredListFileOutputRequiredCompanionsDictionary = TestWDLTool.TestWDLTool_requiredListFileOutputRequiredCompanionsDictionary
    Array[File] TestWDLToolrequiredListFileOutputRequiredCompanionsIndex = TestWDLTool.TestWDLTool_requiredListFileOutputRequiredCompanionsIndex
    Array[File] TestWDLToolrequiredListFileOutputRequiredCompanionsDictionary = TestWDLTool.TestWDLTool_requiredListFileOutputRequiredCompanionsDictionary
    Array[File] TestWDLToolrequiredListFileOutputRequiredCompanionsIndex = TestWDLTool.TestWDLTool_requiredListFileOutputRequiredCompanionsIndex
    Array[File] TestWDLToolrequiredListFileOutputOptionalCompanions = TestWDLTool.TestWDLTool_requiredListFileOutputOptionalCompanions
    Array[File]? TestWDLToolrequiredListFileOutputOptionalCompanionsDictionary = TestWDLTool.TestWDLTool_requiredListFileOutputOptionalCompanionsDictionary
    Array[File]? TestWDLToolrequiredListFileOutputOptionalCompanionsIndex = TestWDLTool.TestWDLTool_requiredListFileOutputOptionalCompanionsIndex
    Array[File] TestWDLToolrequiredListFileOutputMixedCompanions = TestWDLTool.TestWDLTool_requiredListFileOutputMixedCompanions
    Array[File] TestWDLToolrequiredListFileOutputMixedCompanionsRequired = TestWDLTool.TestWDLTool_requiredListFileOutputMixedCompanionsRequired
    Array[File]? TestWDLToolrequiredListFileOutputMixedCompanionsOptional = TestWDLTool.TestWDLTool_requiredListFileOutputMixedCompanionsOptional
    Array[File] TestWDLToolrequiredListFileOutputMixedCompanionsRequired = TestWDLTool.TestWDLTool_requiredListFileOutputMixedCompanionsRequired
  }

  parameter_meta {
    dockerImage: { description: "Docker image for this task" }
    gatk: { description: "Location of gatk to run for this task" }
    memoryRequirements: { description: "Runtime memory requirements for this task" }
    diskRequirements: { description: "Runtime disk requirements for this task" }
    cpuRequirements: { description: "Runtime CPU count for this task" }
    preemptibleRequirements: { description: "Runtime preemptible count for this task" }
    bootdisksizegbRequirements: { description: "Runtime boot disk size for this task" }

    # Positional Arguments
    positionalArgs: { description: "Positional args doc" }
    posDictionary: { description: "Companion resource for positionalArgs" }
    posIndex: { description: "Companion resource for positionalArgs" }

    # Required Arguments
    requiredListFileInputMixedCompanions: { description: "requiredListFileInputMixedCompanions doc" }
    requiredListFileInputMixedCompanionsRequired: { description: "Companion resource for requiredListFileInputMixedCompanions" }
    requiredListFileInputMixedCompanionsOptional: { description: "Companion resource for requiredListFileInputMixedCompanions" }
    requiredListFileInputNoCompanions: { description: "requiredListFileInputNoCompanions doc" }
    requiredListFileInputOptionalCompanions: {
        description: "requiredListFileInputOptionalCompanions doc",
        localization_optional : true 
    }
    requiredListFileInputOptionalCompanionsDictionary: {
      description: "Companion resource for requiredListFileInputOptionalCompanions",
      localization_optional : true 
    }
    requiredListFileInputOptionalCompanionsIndex: {
      description: "Companion resource for requiredListFileInputOptionalCompanions",
      localization_optional : true 
    }
    requiredListFileInputRequiredCompanions: {
        description: "requiredListFileInputRequiredCompanions doc",
        localization_optional : true 
    }
    requiredListFileInputRequiredCompanionsDictionary: {
      description: "Companion resource for requiredListFileInputRequiredCompanions",
      localization_optional : true 
    }
    requiredListFileInputRequiredCompanionsIndex: {
      description: "Companion resource for requiredListFileInputRequiredCompanions",
      localization_optional : true 
    }
    requiredListFileOutputMixedCompanions: { description: "requiredListFileOutputMixedCompanions doc" }
    requiredListFileOutputMixedCompanionsRequired: { description: "Companion resource for requiredListFileOutputMixedCompanions" }
    requiredListFileOutputMixedCompanionsOptional: { description: "Companion resource for requiredListFileOutputMixedCompanions" }
    requiredListFileOutputNoCompanions: { description: "requiredListFileOutputNoCompanions doc" }
    requiredListFileOutputOptionalCompanions: { description: "requiredListFileOutputOptionalCompanions doc" }
    requiredListFileOutputOptionalCompanionsDictionary: { description: "Companion resource for requiredListFileOutputOptionalCompanions" }
    requiredListFileOutputOptionalCompanionsIndex: { description: "Companion resource for requiredListFileOutputOptionalCompanions" }
    requiredListFileOutputRequiredCompanions: { description: "requiredListFileOutputRequiredCompanions doc" }
    requiredListFileOutputRequiredCompanionsDictionary: { description: "Companion resource for requiredListFileOutputRequiredCompanions" }
    requiredListFileOutputRequiredCompanionsIndex: { description: "Companion resource for requiredListFileOutputRequiredCompanions" }
    requiredScalarFileInputNoCompanions: { description: "requiredScalarFileInputNoCompanions doc" }
    requiredScalarFileInputOptionalCompanions: {
        description: "requiredScalarFileInputOptionalCompanions doc",
        localization_optional : true 
    }
    requiredScalarFileInputOptionalCompanionsDictionary: {
      description: "Companion resource for requiredScalarFileInputOptionalCompanions",
      localization_optional : true 
    }
    requiredScalarFileInputOptionalCompanionsIndex: {
      description: "Companion resource for requiredScalarFileInputOptionalCompanions",
      localization_optional : true 
    }
    requiredScalarFileInputRequiredCompanions: {
        description: "requiredScalarFileInputRequiredCompanions doc",
        localization_optional : true 
    }
    requiredScalarFileInputRequiredCompanionsDictionary: {
      description: "Companion resource for requiredScalarFileInputRequiredCompanions",
      localization_optional : true 
    }
    requiredScalarFileInputRequiredCompanionsIndex: {
      description: "Companion resource for requiredScalarFileInputRequiredCompanions",
      localization_optional : true 
    }
    requiredScalarFileOutputNoCompanions: { description: "requiredScalarFileOutputNoCompanions doc" }
    requiredScalarFileOutputOptionalCompanions: { description: "requiredScalarFileOutputOptionalCompanions doc" }
    requiredScalarFileOutputOptionalCompanionsDictionary: { description: "Companion resource for requiredScalarFileOutputOptionalCompanions" }
    requiredScalarFileOutputOptionalCompanionsIndex: { description: "Companion resource for requiredScalarFileOutputOptionalCompanions" }
    requiredScalarFileOutputRequiredCompanions: { description: "requiredScalarFileOutputRequiredCompanions doc" }
    requiredScalarFileOutputRequiredCompanionsDictionary: { description: "Companion resource for requiredScalarFileOutputRequiredCompanions" }
    requiredScalarFileOutputRequiredCompanionsIndex: { description: "Companion resource for requiredScalarFileOutputRequiredCompanions" }
  }
}

task TestWDLTool {

  input {
    String dockerImage
    String gatk
    String memoryRequirements
    String diskRequirements
    String cpuRequirements
    String preemptibleRequirements
    String bootdisksizegbRequirements
    Array[File] positionalArgs
    Array[File] requiredListFileInputMixedCompanions
    Array[File] requiredListFileInputMixedCompanionsRequired
    Array[File]? requiredListFileInputMixedCompanionsOptional
    Array[File] requiredListFileInputNoCompanions
    Array[File] requiredListFileInputOptionalCompanions
    Array[File]? requiredListFileInputOptionalCompanionsDictionary
    Array[File]? requiredListFileInputOptionalCompanionsIndex
    Array[File] requiredListFileInputRequiredCompanions
    Array[File] requiredListFileInputRequiredCompanionsDictionary
    Array[File] requiredListFileInputRequiredCompanionsIndex
    Array[String] requiredListFileOutputMixedCompanions
    Array[String] requiredListFileOutputMixedCompanionsRequired
    Array[String]? requiredListFileOutputMixedCompanionsOptional
    Array[String] requiredListFileOutputNoCompanions
    Array[String] requiredListFileOutputOptionalCompanions
    Array[String]? requiredListFileOutputOptionalCompanionsDictionary
    Array[String]? requiredListFileOutputOptionalCompanionsIndex
    Array[String] requiredListFileOutputRequiredCompanions
    Array[String] requiredListFileOutputRequiredCompanionsDictionary
    Array[String] requiredListFileOutputRequiredCompanionsIndex
    File requiredScalarFileInputNoCompanions
    File requiredScalarFileInputOptionalCompanions
    File? requiredScalarFileInputOptionalCompanionsDictionary
    File? requiredScalarFileInputOptionalCompanionsIndex
    File requiredScalarFileInputRequiredCompanions
    File requiredScalarFileInputRequiredCompanionsDictionary
    File requiredScalarFileInputRequiredCompanionsIndex
    String requiredScalarFileOutputNoCompanions
    String requiredScalarFileOutputOptionalCompanions
    String? requiredScalarFileOutputOptionalCompanionsDictionary
    String? requiredScalarFileOutputOptionalCompanionsIndex
    String requiredScalarFileOutputRequiredCompanions
    String requiredScalarFileOutputRequiredCompanionsDictionary
    String requiredScalarFileOutputRequiredCompanionsIndex

  }

  command <<<
    ~{gatk} TestWDLTool \
    ~{sep=' ' positionalArgs} \
    --requiredListFileInputMixedCompanions ~{sep=' --requiredListFileInputMixedCompanions ' requiredListFileInputMixedCompanions} \
    --requiredListFileInputNoCompanions ~{sep=' --requiredListFileInputNoCompanions ' requiredListFileInputNoCompanions} \
    --requiredListFileInputOptionalCompanions ~{sep=' --requiredListFileInputOptionalCompanions ' requiredListFileInputOptionalCompanions} \
    --requiredListFileInputRequiredCompanions ~{sep=' --requiredListFileInputRequiredCompanions ' requiredListFileInputRequiredCompanions} \
    --requiredListFileOutputMixedCompanions ~{sep=' --requiredListFileOutputMixedCompanions ' requiredListFileOutputMixedCompanions} \
    --requiredListFileOutputNoCompanions ~{sep=' --requiredListFileOutputNoCompanions ' requiredListFileOutputNoCompanions} \
    --requiredListFileOutputOptionalCompanions ~{sep=' --requiredListFileOutputOptionalCompanions ' requiredListFileOutputOptionalCompanions} \
    --requiredListFileOutputRequiredCompanions ~{sep=' --requiredListFileOutputRequiredCompanions ' requiredListFileOutputRequiredCompanions} \
    --requiredScalarFileInputNoCompanions ~{sep=' --requiredScalarFileInputNoCompanions ' requiredScalarFileInputNoCompanions} \
    --requiredScalarFileInputOptionalCompanions ~{sep=' --requiredScalarFileInputOptionalCompanions ' requiredScalarFileInputOptionalCompanions} \
    --requiredScalarFileInputRequiredCompanions ~{sep=' --requiredScalarFileInputRequiredCompanions ' requiredScalarFileInputRequiredCompanions} \
    --requiredScalarFileOutputNoCompanions ~{sep=' --requiredScalarFileOutputNoCompanions ' requiredScalarFileOutputNoCompanions} \
    --requiredScalarFileOutputOptionalCompanions ~{sep=' --requiredScalarFileOutputOptionalCompanions ' requiredScalarFileOutputOptionalCompanions} \
    --requiredScalarFileOutputRequiredCompanions ~{sep=' --requiredScalarFileOutputRequiredCompanions ' requiredScalarFileOutputRequiredCompanions} \

  >>>

  runtime {
      docker: dockerImage
      memory: memoryRequirements
      disks: diskRequirements
      cpu: cpuRequirements
      preemptible: preemptibleRequirements
      bootDiskSizeGb: bootdisksizegbRequirements
  }

  output {
    # Task Outputs                                      
    File TestWDLTool_requiredScalarFileOutputNoCompanions = requiredScalarFileOutputNoCompanions
    File TestWDLTool_requiredScalarFileOutputRequiredCompanions = requiredScalarFileOutputRequiredCompanions
    File TestWDLTool_requiredScalarFileOutputRequiredCompanionsDictionary = requiredScalarFileOutputRequiredCompanionsDictionary
    File TestWDLTool_requiredScalarFileOutputRequiredCompanionsIndex = requiredScalarFileOutputRequiredCompanionsIndex
    File TestWDLTool_requiredScalarFileOutputRequiredCompanionsDictionary = requiredScalarFileOutputRequiredCompanionsDictionary
    File TestWDLTool_requiredScalarFileOutputRequiredCompanionsIndex = requiredScalarFileOutputRequiredCompanionsIndex
    File TestWDLTool_requiredScalarFileOutputOptionalCompanions = requiredScalarFileOutputOptionalCompanions
    File? TestWDLTool_requiredScalarFileOutputOptionalCompanionsDictionary = requiredScalarFileOutputOptionalCompanionsDictionary
    File? TestWDLTool_requiredScalarFileOutputOptionalCompanionsIndex = requiredScalarFileOutputOptionalCompanionsIndex
    Array[File] TestWDLTool_requiredListFileOutputNoCompanions = requiredListFileOutputNoCompanions
    Array[File] TestWDLTool_requiredListFileOutputRequiredCompanions = requiredListFileOutputRequiredCompanions
    Array[File] TestWDLTool_requiredListFileOutputRequiredCompanionsDictionary = requiredListFileOutputRequiredCompanionsDictionary
    Array[File] TestWDLTool_requiredListFileOutputRequiredCompanionsIndex = requiredListFileOutputRequiredCompanionsIndex
    Array[File] TestWDLTool_requiredListFileOutputRequiredCompanionsDictionary = requiredListFileOutputRequiredCompanionsDictionary
    Array[File] TestWDLTool_requiredListFileOutputRequiredCompanionsIndex = requiredListFileOutputRequiredCompanionsIndex
    Array[File] TestWDLTool_requiredListFileOutputOptionalCompanions = requiredListFileOutputOptionalCompanions
    Array[File]? TestWDLTool_requiredListFileOutputOptionalCompanionsDictionary = requiredListFileOutputOptionalCompanionsDictionary
    Array[File]? TestWDLTool_requiredListFileOutputOptionalCompanionsIndex = requiredListFileOutputOptionalCompanionsIndex
    Array[File] TestWDLTool_requiredListFileOutputMixedCompanions = requiredListFileOutputMixedCompanions
    Array[File] TestWDLTool_requiredListFileOutputMixedCompanionsRequired = requiredListFileOutputMixedCompanionsRequired
    Array[File]? TestWDLTool_requiredListFileOutputMixedCompanionsOptional = requiredListFileOutputMixedCompanionsOptional
    Array[File] TestWDLTool_requiredListFileOutputMixedCompanionsRequired = requiredListFileOutputMixedCompanionsRequired
  }

  parameter_meta {
    dockerImage: { description: "Docker image for this task" }
    gatk: { description: "Location of gatk to run for this task" }
    memoryRequirements: { description: "Runtime memory requirements for this task" }
    diskRequirements: { description: "Runtime disk requirements for this task" }
    cpuRequirements: { description: "Runtime CPU count for this task" }
    preemptibleRequirements: { description: "Runtime preemptible count for this task" }
    bootdisksizegbRequirements: { description: "Runtime boot disk size for this task" }

    # Positional Arguments
    positionalArgs: { description: "Positional args doc" }
    posDictionary: { description: "Companion resource for positionalArgs" }
    posIndex: { description: "Companion resource for positionalArgs" }

    # Required Arguments
    requiredListFileInputMixedCompanions: { description: "requiredListFileInputMixedCompanions doc" }
    requiredListFileInputMixedCompanionsRequired: { description: "Companion resource for requiredListFileInputMixedCompanions" }
    requiredListFileInputMixedCompanionsOptional: { description: "Companion resource for requiredListFileInputMixedCompanions" }
    requiredListFileInputNoCompanions: { description: "requiredListFileInputNoCompanions doc" }
    requiredListFileInputOptionalCompanions: {
        description: "requiredListFileInputOptionalCompanions doc",
        localization_optional : true 
    }
    requiredListFileInputOptionalCompanionsDictionary: {
      description: "Companion resource for requiredListFileInputOptionalCompanions",
      localization_optional : true 
    }
    requiredListFileInputOptionalCompanionsIndex: {
      description: "Companion resource for requiredListFileInputOptionalCompanions",
      localization_optional : true 
    }
    requiredListFileInputRequiredCompanions: {
        description: "requiredListFileInputRequiredCompanions doc",
        localization_optional : true 
    }
    requiredListFileInputRequiredCompanionsDictionary: {
      description: "Companion resource for requiredListFileInputRequiredCompanions",
      localization_optional : true 
    }
    requiredListFileInputRequiredCompanionsIndex: {
      description: "Companion resource for requiredListFileInputRequiredCompanions",
      localization_optional : true 
    }
    requiredListFileOutputMixedCompanions: { description: "requiredListFileOutputMixedCompanions doc" }
    requiredListFileOutputMixedCompanionsRequired: { description: "Companion resource for requiredListFileOutputMixedCompanions" }
    requiredListFileOutputMixedCompanionsOptional: { description: "Companion resource for requiredListFileOutputMixedCompanions" }
    requiredListFileOutputNoCompanions: { description: "requiredListFileOutputNoCompanions doc" }
    requiredListFileOutputOptionalCompanions: { description: "requiredListFileOutputOptionalCompanions doc" }
    requiredListFileOutputOptionalCompanionsDictionary: { description: "Companion resource for requiredListFileOutputOptionalCompanions" }
    requiredListFileOutputOptionalCompanionsIndex: { description: "Companion resource for requiredListFileOutputOptionalCompanions" }
    requiredListFileOutputRequiredCompanions: { description: "requiredListFileOutputRequiredCompanions doc" }
    requiredListFileOutputRequiredCompanionsDictionary: { description: "Companion resource for requiredListFileOutputRequiredCompanions" }
    requiredListFileOutputRequiredCompanionsIndex: { description: "Companion resource for requiredListFileOutputRequiredCompanions" }
    requiredScalarFileInputNoCompanions: { description: "requiredScalarFileInputNoCompanions doc" }
    requiredScalarFileInputOptionalCompanions: {
        description: "requiredScalarFileInputOptionalCompanions doc",
        localization_optional : true 
    }
    requiredScalarFileInputOptionalCompanionsDictionary: {
      description: "Companion resource for requiredScalarFileInputOptionalCompanions",
      localization_optional : true 
    }
    requiredScalarFileInputOptionalCompanionsIndex: {
      description: "Companion resource for requiredScalarFileInputOptionalCompanions",
      localization_optional : true 
    }
    requiredScalarFileInputRequiredCompanions: {
        description: "requiredScalarFileInputRequiredCompanions doc",
        localization_optional : true 
    }
    requiredScalarFileInputRequiredCompanionsDictionary: {
      description: "Companion resource for requiredScalarFileInputRequiredCompanions",
      localization_optional : true 
    }
    requiredScalarFileInputRequiredCompanionsIndex: {
      description: "Companion resource for requiredScalarFileInputRequiredCompanions",
      localization_optional : true 
    }
    requiredScalarFileOutputNoCompanions: { description: "requiredScalarFileOutputNoCompanions doc" }
    requiredScalarFileOutputOptionalCompanions: { description: "requiredScalarFileOutputOptionalCompanions doc" }
    requiredScalarFileOutputOptionalCompanionsDictionary: { description: "Companion resource for requiredScalarFileOutputOptionalCompanions" }
    requiredScalarFileOutputOptionalCompanionsIndex: { description: "Companion resource for requiredScalarFileOutputOptionalCompanions" }
    requiredScalarFileOutputRequiredCompanions: { description: "requiredScalarFileOutputRequiredCompanions doc" }
    requiredScalarFileOutputRequiredCompanionsDictionary: { description: "Companion resource for requiredScalarFileOutputRequiredCompanions" }
    requiredScalarFileOutputRequiredCompanionsIndex: { description: "Companion resource for requiredScalarFileOutputRequiredCompanions" }
  }
}

