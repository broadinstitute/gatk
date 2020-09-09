version 1.0

# Run TestWDLTool (WDL auto generated from: GATK Version 1.2.3.4)
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
#    requiredListFileInputNoCompanions                  requiredListFileInputNoCompanions doc                       
#    requiredListFileInputOptionalCompanions            requiredListFileInputOptionalCompanions doc                 
#    requiredListFileInputOptionalCompanionsDictionary  Optional Companion resource for requiredListFileInputOptionalCompanions
#    requiredListFileInputOptionalCompanionsIndex       Optional Companion resource for requiredListFileInputOptionalCompanions
#    requiredListFileInputRequiredCompanions            requiredListFileInputRequiredCompanions doc                 
#    requiredListFileInputRequiredCompanionsDictionary  Companion resource for requiredListFileInputRequiredCompanions
#    requiredListFileInputRequiredCompanionsIndex       Companion resource for requiredListFileInputRequiredCompanions
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
#  Optional Tool Arguments
#    optionalListDoubleInput                            optionalListDoubleInput doc                                 
#    optionalListFileInputNoCompanions                  optionalListFileInputNoCompanions doc                       
#    optionalListFileInputOptionalCompanions            optionalListFileInputOptionalCompanions doc                 
#    optionalListFileInputOptionalCompanionsDictionary  Optional Companion resource for optionalListFileInputOptionalCompanions
#    optionalListFileInputOptionalCompanionsIndex       Optional Companion resource for optionalListFileInputOptionalCompanions
#    optionalListFileInputRequiredCompanions            optionalListFileInputRequiredCompanions doc                 
#    optionalListFileInputRequiredCompanionsDictionary  Optional Companion resource for optionalListFileInputRequiredCompanions
#    optionalListFileInputRequiredCompanionsIndex       Optional Companion resource for optionalListFileInputRequiredCompanions
#    optionalListFileOutputRequiredCompanions           optionalListFileOutputRequiredCompanions doc                
#    optionalListFileOutputRequiredCompanionsDictionary Optional Companion resource for optionalListFileOutputRequiredCompanions
#    optionalListFileOutputRequiredCompanionsIndex      Optional Companion resource for optionalListFileOutputRequiredCompanions
#    optionalListFloatInput                             optionalListFloatInput doc                                  
#    optionalListIntegerInput                           optionalListIntegerInput doc                                
#    optionalListLongInput                              optionalListLongInput doc                                   
#    optionalListStringInput                            optionalListStringInput doc                                 
#    optionalScalarDoubleInput                          optionalScalarDoubleInput doc                               
#    optionalScalarDoublePrimitiveInput                 optionalScalarDoublePrimitiveInput doc                      
#    optionalScalarFileInputNoCompanions                optionalScalarFileInputNoCompanions doc                     
#    optionalScalarFileInputOptionalCompanions          optionalScalarFileInputOptionalCompanions doc               
#    optionalScalarFileInputOptionalCompanionsDictionary Optional Companion resource for optionalScalarFileInputOptionalCompanions
#    optionalScalarFileInputOptionalCompanionsIndex     Optional Companion resource for optionalScalarFileInputOptionalCompanions
#    optionalScalarFileInputRequiredCompanions          optionalScalarFileInputRequiredCompanions doc               
#    optionalScalarFileInputRequiredCompanionsDictionary Optional Companion resource for optionalScalarFileInputRequiredCompanions
#    optionalScalarFileInputRequiredCompanionsIndex     Optional Companion resource for optionalScalarFileInputRequiredCompanions
#    optionalScalarFileOutputNoCompanions               optionalScalarFileOutputNoCompanions doc                    
#    optionalScalarFileOutputOptionalCompanions         optionalScalarFileOutputOptionalCompanions doc              
#    optionalScalarFileOutputOptionalCompanionsDictionary Optional Companion resource for optionalScalarFileOutputOptionalCompanions
#    optionalScalarFileOutputOptionalCompanionsIndex    Optional Companion resource for optionalScalarFileOutputOptionalCompanions
#    optionalScalarFileOutputRequiredCompanions         optionalScalarFileOutputRequiredCompanions doc              
#    optionalScalarFileOutputRequiredCompanionsDictionary Optional Companion resource for optionalScalarFileOutputRequiredCompanions
#    optionalScalarFileOutputRequiredCompanionsIndex    Optional Companion resource for optionalScalarFileOutputRequiredCompanions
#    optionalScalarFloatInput                           optionalScalarFloatInput doc                                
#    optionalScalarFloatPrimitiveInput                  optionalScalarFloatPrimitiveInput doc                       
#    optionalScalarIntegerInput                         optionalScalarIntegerInput doc                              
#    optionalScalarIntegerPrimitiveInput                optionalScalarIntegerPrimitiveInput doc                     
#    optionalScalarLongInput                            optionalScalarLongInput doc                                 
#    optionalScalarLongPrimitiveInput                   optionalScalarLongPrimitiveInput doc                        
#    optionalScalarStringInput                          optionalScalarStringInput doc                               
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
    Array[File] requiredListFileInputNoCompanions
    Array[File] requiredListFileInputOptionalCompanions
    Array[File]? requiredListFileInputOptionalCompanionsDictionary
    Array[File]? requiredListFileInputOptionalCompanionsIndex
    Array[File] requiredListFileInputRequiredCompanions
    Array[File] requiredListFileInputRequiredCompanionsDictionary
    Array[File] requiredListFileInputRequiredCompanionsIndex
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

    # Optional Tool Arguments
    Array[Float]? optionalListDoubleInput
    Array[File]? optionalListFileInputNoCompanions
    Array[File]? optionalListFileInputOptionalCompanions
    Array[File]? optionalListFileInputOptionalCompanionsDictionary
    Array[File]? optionalListFileInputOptionalCompanionsIndex
    Array[File]? optionalListFileInputRequiredCompanions
    Array[File]? optionalListFileInputRequiredCompanionsDictionary
    Array[File]? optionalListFileInputRequiredCompanionsIndex
    Array[String]? optionalListFileOutputRequiredCompanions
    Array[String]? optionalListFileOutputRequiredCompanionsDictionary
    Array[String]? optionalListFileOutputRequiredCompanionsIndex
    Array[Float]? optionalListFloatInput
    Array[Int]? optionalListIntegerInput
    Array[Int]? optionalListLongInput
    Array[String]? optionalListStringInput
    Float? optionalScalarDoubleInput
    Float? optionalScalarDoublePrimitiveInput
    File? optionalScalarFileInputNoCompanions
    File? optionalScalarFileInputOptionalCompanions
    File? optionalScalarFileInputOptionalCompanionsDictionary
    File? optionalScalarFileInputOptionalCompanionsIndex
    File? optionalScalarFileInputRequiredCompanions
    File? optionalScalarFileInputRequiredCompanionsDictionary
    File? optionalScalarFileInputRequiredCompanionsIndex
    String? optionalScalarFileOutputNoCompanions
    String? optionalScalarFileOutputOptionalCompanions
    String? optionalScalarFileOutputOptionalCompanionsDictionary
    String? optionalScalarFileOutputOptionalCompanionsIndex
    String? optionalScalarFileOutputRequiredCompanions
    String? optionalScalarFileOutputRequiredCompanionsDictionary
    String? optionalScalarFileOutputRequiredCompanionsIndex
    Float? optionalScalarFloatInput
    Float? optionalScalarFloatPrimitiveInput
    Int? optionalScalarIntegerInput
    Int? optionalScalarIntegerPrimitiveInput
    Int? optionalScalarLongInput
    Int? optionalScalarLongPrimitiveInput
    String? optionalScalarStringInput

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
        requiredListFileInputNoCompanions                  = requiredListFileInputNoCompanions,
        requiredListFileInputOptionalCompanions            = requiredListFileInputOptionalCompanions,
        requiredListFileInputOptionalCompanionsDictionary  = requiredListFileInputOptionalCompanionsDictionary,
        requiredListFileInputOptionalCompanionsIndex       = requiredListFileInputOptionalCompanionsIndex,
        requiredListFileInputRequiredCompanions            = requiredListFileInputRequiredCompanions,
        requiredListFileInputRequiredCompanionsDictionary  = requiredListFileInputRequiredCompanionsDictionary,
        requiredListFileInputRequiredCompanionsIndex       = requiredListFileInputRequiredCompanionsIndex,
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

        # Optional Tool Arguments
        optionalListDoubleInput                            = optionalListDoubleInput,
        optionalListFileInputNoCompanions                  = optionalListFileInputNoCompanions,
        optionalListFileInputOptionalCompanions            = optionalListFileInputOptionalCompanions,
        optionalListFileInputOptionalCompanionsDictionary  = optionalListFileInputOptionalCompanionsDictionary,
        optionalListFileInputOptionalCompanionsIndex       = optionalListFileInputOptionalCompanionsIndex,
        optionalListFileInputRequiredCompanions            = optionalListFileInputRequiredCompanions,
        optionalListFileInputRequiredCompanionsDictionary  = optionalListFileInputRequiredCompanionsDictionary,
        optionalListFileInputRequiredCompanionsIndex       = optionalListFileInputRequiredCompanionsIndex,
        optionalListFileOutputRequiredCompanions           = optionalListFileOutputRequiredCompanions,
        optionalListFileOutputRequiredCompanionsDictionary = optionalListFileOutputRequiredCompanionsDictionary,
        optionalListFileOutputRequiredCompanionsIndex      = optionalListFileOutputRequiredCompanionsIndex,
        optionalListFloatInput                             = optionalListFloatInput,
        optionalListIntegerInput                           = optionalListIntegerInput,
        optionalListLongInput                              = optionalListLongInput,
        optionalListStringInput                            = optionalListStringInput,
        optionalScalarDoubleInput                          = optionalScalarDoubleInput,
        optionalScalarDoublePrimitiveInput                 = optionalScalarDoublePrimitiveInput,
        optionalScalarFileInputNoCompanions                = optionalScalarFileInputNoCompanions,
        optionalScalarFileInputOptionalCompanions          = optionalScalarFileInputOptionalCompanions,
        optionalScalarFileInputOptionalCompanionsDictionary = optionalScalarFileInputOptionalCompanionsDictionary,
        optionalScalarFileInputOptionalCompanionsIndex     = optionalScalarFileInputOptionalCompanionsIndex,
        optionalScalarFileInputRequiredCompanions          = optionalScalarFileInputRequiredCompanions,
        optionalScalarFileInputRequiredCompanionsDictionary = optionalScalarFileInputRequiredCompanionsDictionary,
        optionalScalarFileInputRequiredCompanionsIndex     = optionalScalarFileInputRequiredCompanionsIndex,
        optionalScalarFileOutputNoCompanions               = optionalScalarFileOutputNoCompanions,
        optionalScalarFileOutputOptionalCompanions         = optionalScalarFileOutputOptionalCompanions,
        optionalScalarFileOutputOptionalCompanionsDictionary = optionalScalarFileOutputOptionalCompanionsDictionary,
        optionalScalarFileOutputOptionalCompanionsIndex    = optionalScalarFileOutputOptionalCompanionsIndex,
        optionalScalarFileOutputRequiredCompanions         = optionalScalarFileOutputRequiredCompanions,
        optionalScalarFileOutputRequiredCompanionsDictionary = optionalScalarFileOutputRequiredCompanionsDictionary,
        optionalScalarFileOutputRequiredCompanionsIndex    = optionalScalarFileOutputRequiredCompanionsIndex,
        optionalScalarFloatInput                           = optionalScalarFloatInput,
        optionalScalarFloatPrimitiveInput                  = optionalScalarFloatPrimitiveInput,
        optionalScalarIntegerInput                         = optionalScalarIntegerInput,
        optionalScalarIntegerPrimitiveInput                = optionalScalarIntegerPrimitiveInput,
        optionalScalarLongInput                            = optionalScalarLongInput,
        optionalScalarLongPrimitiveInput                   = optionalScalarLongPrimitiveInput,
        optionalScalarStringInput                          = optionalScalarStringInput,

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

    # Optional Tool Arguments
    optionalListDoubleInput: { description: "optionalListDoubleInput doc" }
    optionalListFileInputNoCompanions: { description: "optionalListFileInputNoCompanions doc" }
    optionalListFileInputOptionalCompanions: { description: "optionalListFileInputOptionalCompanions doc" }
    optionalListFileInputOptionalCompanionsDictionary: { description: "Companion resource for optionalListFileInputOptionalCompanions" }
    optionalListFileInputOptionalCompanionsIndex: { description: "Companion resource for optionalListFileInputOptionalCompanions" }
    optionalListFileInputRequiredCompanions: { description: "optionalListFileInputRequiredCompanions doc" }
    optionalListFileInputRequiredCompanionsDictionary: { description: "Companion resource for optionalListFileInputRequiredCompanions" }
    optionalListFileInputRequiredCompanionsIndex: { description: "Companion resource for optionalListFileInputRequiredCompanions" }
    optionalListFileOutputRequiredCompanions: { description: "optionalListFileOutputRequiredCompanions doc" }
    optionalListFileOutputRequiredCompanionsDictionary: { description: "Companion resource for optionalListFileOutputRequiredCompanions" }
    optionalListFileOutputRequiredCompanionsIndex: { description: "Companion resource for optionalListFileOutputRequiredCompanions" }
    optionalListFloatInput: { description: "optionalListFloatInput doc" }
    optionalListIntegerInput: { description: "optionalListIntegerInput doc" }
    optionalListLongInput: { description: "optionalListLongInput doc" }
    optionalListStringInput: { description: "optionalListStringInput doc" }
    optionalScalarDoubleInput: { description: "optionalScalarDoubleInput doc" }
    optionalScalarDoublePrimitiveInput: { description: "optionalScalarDoublePrimitiveInput doc" }
    optionalScalarFileInputNoCompanions: { description: "optionalScalarFileInputNoCompanions doc" }
    optionalScalarFileInputOptionalCompanions: { description: "optionalScalarFileInputOptionalCompanions doc" }
    optionalScalarFileInputOptionalCompanionsDictionary: { description: "Companion resource for optionalScalarFileInputOptionalCompanions" }
    optionalScalarFileInputOptionalCompanionsIndex: { description: "Companion resource for optionalScalarFileInputOptionalCompanions" }
    optionalScalarFileInputRequiredCompanions: { description: "optionalScalarFileInputRequiredCompanions doc" }
    optionalScalarFileInputRequiredCompanionsDictionary: { description: "Companion resource for optionalScalarFileInputRequiredCompanions" }
    optionalScalarFileInputRequiredCompanionsIndex: { description: "Companion resource for optionalScalarFileInputRequiredCompanions" }
    optionalScalarFileOutputNoCompanions: { description: "optionalScalarFileOutputNoCompanions doc" }
    optionalScalarFileOutputOptionalCompanions: { description: "optionalScalarFileOutputOptionalCompanions doc" }
    optionalScalarFileOutputOptionalCompanionsDictionary: { description: "Companion resource for optionalScalarFileOutputOptionalCompanions" }
    optionalScalarFileOutputOptionalCompanionsIndex: { description: "Companion resource for optionalScalarFileOutputOptionalCompanions" }
    optionalScalarFileOutputRequiredCompanions: { description: "optionalScalarFileOutputRequiredCompanions doc" }
    optionalScalarFileOutputRequiredCompanionsDictionary: { description: "Companion resource for optionalScalarFileOutputRequiredCompanions" }
    optionalScalarFileOutputRequiredCompanionsIndex: { description: "Companion resource for optionalScalarFileOutputRequiredCompanions" }
    optionalScalarFloatInput: { description: "optionalScalarFloatInput doc" }
    optionalScalarFloatPrimitiveInput: { description: "optionalScalarFloatPrimitiveInput doc" }
    optionalScalarIntegerInput: { description: "optionalScalarIntegerInput doc" }
    optionalScalarIntegerPrimitiveInput: { description: "optionalScalarIntegerPrimitiveInput doc" }
    optionalScalarLongInput: { description: "optionalScalarLongInput doc" }
    optionalScalarLongPrimitiveInput: { description: "optionalScalarLongPrimitiveInput doc" }
    optionalScalarStringInput: { description: "optionalScalarStringInput doc" }
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
    Array[File] requiredListFileInputNoCompanions
    Array[File] requiredListFileInputOptionalCompanions
    Array[File]? requiredListFileInputOptionalCompanionsDictionary
    Array[File]? requiredListFileInputOptionalCompanionsIndex
    Array[File] requiredListFileInputRequiredCompanions
    Array[File] requiredListFileInputRequiredCompanionsDictionary
    Array[File] requiredListFileInputRequiredCompanionsIndex
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
    Array[Float]? optionalListDoubleInput
    Array[File]? optionalListFileInputNoCompanions
    Array[File]? optionalListFileInputOptionalCompanions
    Array[File]? optionalListFileInputOptionalCompanionsDictionary
    Array[File]? optionalListFileInputOptionalCompanionsIndex
    Array[File]? optionalListFileInputRequiredCompanions
    Array[File]? optionalListFileInputRequiredCompanionsDictionary
    Array[File]? optionalListFileInputRequiredCompanionsIndex
    Array[String]? optionalListFileOutputRequiredCompanions
    Array[String]? optionalListFileOutputRequiredCompanionsDictionary
    Array[String]? optionalListFileOutputRequiredCompanionsIndex
    Array[Float]? optionalListFloatInput
    Array[Int]? optionalListIntegerInput
    Array[Int]? optionalListLongInput
    Array[String]? optionalListStringInput
    Float? optionalScalarDoubleInput
    Float? optionalScalarDoublePrimitiveInput
    File? optionalScalarFileInputNoCompanions
    File? optionalScalarFileInputOptionalCompanions
    File? optionalScalarFileInputOptionalCompanionsDictionary
    File? optionalScalarFileInputOptionalCompanionsIndex
    File? optionalScalarFileInputRequiredCompanions
    File? optionalScalarFileInputRequiredCompanionsDictionary
    File? optionalScalarFileInputRequiredCompanionsIndex
    String? optionalScalarFileOutputNoCompanions
    String? optionalScalarFileOutputOptionalCompanions
    String? optionalScalarFileOutputOptionalCompanionsDictionary
    String? optionalScalarFileOutputOptionalCompanionsIndex
    String? optionalScalarFileOutputRequiredCompanions
    String? optionalScalarFileOutputRequiredCompanionsDictionary
    String? optionalScalarFileOutputRequiredCompanionsIndex
    Float? optionalScalarFloatInput
    Float? optionalScalarFloatPrimitiveInput
    Int? optionalScalarIntegerInput
    Int? optionalScalarIntegerPrimitiveInput
    Int? optionalScalarLongInput
    Int? optionalScalarLongPrimitiveInput
    String? optionalScalarStringInput

  }

  command <<<
    ~{gatk} TestWDLTool \
    ~{sep=' ' positionalArgs} \
    --requiredListFileInputNoCompanions ~{sep=' --requiredListFileInputNoCompanions ' requiredListFileInputNoCompanions} \
    --requiredListFileInputOptionalCompanions ~{sep=' --requiredListFileInputOptionalCompanions ' requiredListFileInputOptionalCompanions} \
    --requiredListFileInputRequiredCompanions ~{sep=' --requiredListFileInputRequiredCompanions ' requiredListFileInputRequiredCompanions} \
    --requiredListFileOutputNoCompanions ~{sep=' --requiredListFileOutputNoCompanions ' requiredListFileOutputNoCompanions} \
    --requiredListFileOutputOptionalCompanions ~{sep=' --requiredListFileOutputOptionalCompanions ' requiredListFileOutputOptionalCompanions} \
    --requiredListFileOutputRequiredCompanions ~{sep=' --requiredListFileOutputRequiredCompanions ' requiredListFileOutputRequiredCompanions} \
    --requiredScalarFileInputNoCompanions ~{sep=' --requiredScalarFileInputNoCompanions ' requiredScalarFileInputNoCompanions} \
    --requiredScalarFileInputOptionalCompanions ~{sep=' --requiredScalarFileInputOptionalCompanions ' requiredScalarFileInputOptionalCompanions} \
    --requiredScalarFileInputRequiredCompanions ~{sep=' --requiredScalarFileInputRequiredCompanions ' requiredScalarFileInputRequiredCompanions} \
    --requiredScalarFileOutputNoCompanions ~{sep=' --requiredScalarFileOutputNoCompanions ' requiredScalarFileOutputNoCompanions} \
    --requiredScalarFileOutputOptionalCompanions ~{sep=' --requiredScalarFileOutputOptionalCompanions ' requiredScalarFileOutputOptionalCompanions} \
    --requiredScalarFileOutputRequiredCompanions ~{sep=' --requiredScalarFileOutputRequiredCompanions ' requiredScalarFileOutputRequiredCompanions} \
    ~{true='--optionalListDoubleInput ' false='' defined(optionalListDoubleInput)}~{sep=' --optionalListDoubleInput ' optionalListDoubleInput} \
    ~{true='--optionalListFileInputNoCompanions ' false='' defined(optionalListFileInputNoCompanions)}~{sep=' --optionalListFileInputNoCompanions ' optionalListFileInputNoCompanions} \
    ~{true='--optionalListFileInputOptionalCompanions ' false='' defined(optionalListFileInputOptionalCompanions)}~{sep=' --optionalListFileInputOptionalCompanions ' optionalListFileInputOptionalCompanions} \
    ~{true='--optionalListFileInputRequiredCompanions ' false='' defined(optionalListFileInputRequiredCompanions)}~{sep=' --optionalListFileInputRequiredCompanions ' optionalListFileInputRequiredCompanions} \
    ~{true='--optionalListFileOutputRequiredCompanions ' false='' defined(optionalListFileOutputRequiredCompanions)}~{sep=' --optionalListFileOutputRequiredCompanions ' optionalListFileOutputRequiredCompanions} \
    ~{true='--optionalListFloatInput ' false='' defined(optionalListFloatInput)}~{sep=' --optionalListFloatInput ' optionalListFloatInput} \
    ~{true='--optionalListIntegerInput ' false='' defined(optionalListIntegerInput)}~{sep=' --optionalListIntegerInput ' optionalListIntegerInput} \
    ~{true='--optionalListLongInput ' false='' defined(optionalListLongInput)}~{sep=' --optionalListLongInput ' optionalListLongInput} \
    ~{true='--optionalListStringInput ' false='' defined(optionalListStringInput)}~{sep=' --optionalListStringInput ' optionalListStringInput} \
    ~{true='--optionalScalarDoubleInput ' false='' defined(optionalScalarDoubleInput)}~{sep=' --optionalScalarDoubleInput ' optionalScalarDoubleInput} \
    ~{true='--optionalScalarDoublePrimitiveInput ' false='' defined(optionalScalarDoublePrimitiveInput)}~{sep=' --optionalScalarDoublePrimitiveInput ' optionalScalarDoublePrimitiveInput} \
    ~{true='--optionalScalarFileInputNoCompanions ' false='' defined(optionalScalarFileInputNoCompanions)}~{sep=' --optionalScalarFileInputNoCompanions ' optionalScalarFileInputNoCompanions} \
    ~{true='--optionalScalarFileInputOptionalCompanions ' false='' defined(optionalScalarFileInputOptionalCompanions)}~{sep=' --optionalScalarFileInputOptionalCompanions ' optionalScalarFileInputOptionalCompanions} \
    ~{true='--optionalScalarFileInputRequiredCompanions ' false='' defined(optionalScalarFileInputRequiredCompanions)}~{sep=' --optionalScalarFileInputRequiredCompanions ' optionalScalarFileInputRequiredCompanions} \
    ~{true='--optionalScalarFileOutputNoCompanions ' false='' defined(optionalScalarFileOutputNoCompanions)}~{sep=' --optionalScalarFileOutputNoCompanions ' optionalScalarFileOutputNoCompanions} \
    ~{true='--optionalScalarFileOutputOptionalCompanions ' false='' defined(optionalScalarFileOutputOptionalCompanions)}~{sep=' --optionalScalarFileOutputOptionalCompanions ' optionalScalarFileOutputOptionalCompanions} \
    ~{true='--optionalScalarFileOutputRequiredCompanions ' false='' defined(optionalScalarFileOutputRequiredCompanions)}~{sep=' --optionalScalarFileOutputRequiredCompanions ' optionalScalarFileOutputRequiredCompanions} \
    ~{true='--optionalScalarFloatInput ' false='' defined(optionalScalarFloatInput)}~{sep=' --optionalScalarFloatInput ' optionalScalarFloatInput} \
    ~{true='--optionalScalarFloatPrimitiveInput ' false='' defined(optionalScalarFloatPrimitiveInput)}~{sep=' --optionalScalarFloatPrimitiveInput ' optionalScalarFloatPrimitiveInput} \
    ~{true='--optionalScalarIntegerInput ' false='' defined(optionalScalarIntegerInput)}~{sep=' --optionalScalarIntegerInput ' optionalScalarIntegerInput} \
    ~{true='--optionalScalarIntegerPrimitiveInput ' false='' defined(optionalScalarIntegerPrimitiveInput)}~{sep=' --optionalScalarIntegerPrimitiveInput ' optionalScalarIntegerPrimitiveInput} \
    ~{true='--optionalScalarLongInput ' false='' defined(optionalScalarLongInput)}~{sep=' --optionalScalarLongInput ' optionalScalarLongInput} \
    ~{true='--optionalScalarLongPrimitiveInput ' false='' defined(optionalScalarLongPrimitiveInput)}~{sep=' --optionalScalarLongPrimitiveInput ' optionalScalarLongPrimitiveInput} \
    ~{true='--optionalScalarStringInput ' false='' defined(optionalScalarStringInput)}~{sep=' --optionalScalarStringInput ' optionalScalarStringInput} \

  >>>

  runtime {
      #docker: dockerImage
      memory: memoryRequirements
      disks: diskRequirements
      cpu: cpuRequirements
      preemptible: preemptibleRequirements
      bootDiskSizeGb: bootdisksizegbRequirements
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

    # Optional Tool Arguments
    optionalListDoubleInput: { description: "optionalListDoubleInput doc" }
    optionalListFileInputNoCompanions: { description: "optionalListFileInputNoCompanions doc" }
    optionalListFileInputOptionalCompanions: { description: "optionalListFileInputOptionalCompanions doc" }
    optionalListFileInputOptionalCompanionsDictionary: { description: "Companion resource for optionalListFileInputOptionalCompanions" }
    optionalListFileInputOptionalCompanionsIndex: { description: "Companion resource for optionalListFileInputOptionalCompanions" }
    optionalListFileInputRequiredCompanions: { description: "optionalListFileInputRequiredCompanions doc" }
    optionalListFileInputRequiredCompanionsDictionary: { description: "Companion resource for optionalListFileInputRequiredCompanions" }
    optionalListFileInputRequiredCompanionsIndex: { description: "Companion resource for optionalListFileInputRequiredCompanions" }
    optionalListFileOutputRequiredCompanions: { description: "optionalListFileOutputRequiredCompanions doc" }
    optionalListFileOutputRequiredCompanionsDictionary: { description: "Companion resource for optionalListFileOutputRequiredCompanions" }
    optionalListFileOutputRequiredCompanionsIndex: { description: "Companion resource for optionalListFileOutputRequiredCompanions" }
    optionalListFloatInput: { description: "optionalListFloatInput doc" }
    optionalListIntegerInput: { description: "optionalListIntegerInput doc" }
    optionalListLongInput: { description: "optionalListLongInput doc" }
    optionalListStringInput: { description: "optionalListStringInput doc" }
    optionalScalarDoubleInput: { description: "optionalScalarDoubleInput doc" }
    optionalScalarDoublePrimitiveInput: { description: "optionalScalarDoublePrimitiveInput doc" }
    optionalScalarFileInputNoCompanions: { description: "optionalScalarFileInputNoCompanions doc" }
    optionalScalarFileInputOptionalCompanions: { description: "optionalScalarFileInputOptionalCompanions doc" }
    optionalScalarFileInputOptionalCompanionsDictionary: { description: "Companion resource for optionalScalarFileInputOptionalCompanions" }
    optionalScalarFileInputOptionalCompanionsIndex: { description: "Companion resource for optionalScalarFileInputOptionalCompanions" }
    optionalScalarFileInputRequiredCompanions: { description: "optionalScalarFileInputRequiredCompanions doc" }
    optionalScalarFileInputRequiredCompanionsDictionary: { description: "Companion resource for optionalScalarFileInputRequiredCompanions" }
    optionalScalarFileInputRequiredCompanionsIndex: { description: "Companion resource for optionalScalarFileInputRequiredCompanions" }
    optionalScalarFileOutputNoCompanions: { description: "optionalScalarFileOutputNoCompanions doc" }
    optionalScalarFileOutputOptionalCompanions: { description: "optionalScalarFileOutputOptionalCompanions doc" }
    optionalScalarFileOutputOptionalCompanionsDictionary: { description: "Companion resource for optionalScalarFileOutputOptionalCompanions" }
    optionalScalarFileOutputOptionalCompanionsIndex: { description: "Companion resource for optionalScalarFileOutputOptionalCompanions" }
    optionalScalarFileOutputRequiredCompanions: { description: "optionalScalarFileOutputRequiredCompanions doc" }
    optionalScalarFileOutputRequiredCompanionsDictionary: { description: "Companion resource for optionalScalarFileOutputRequiredCompanions" }
    optionalScalarFileOutputRequiredCompanionsIndex: { description: "Companion resource for optionalScalarFileOutputRequiredCompanions" }
    optionalScalarFloatInput: { description: "optionalScalarFloatInput doc" }
    optionalScalarFloatPrimitiveInput: { description: "optionalScalarFloatPrimitiveInput doc" }
    optionalScalarIntegerInput: { description: "optionalScalarIntegerInput doc" }
    optionalScalarIntegerPrimitiveInput: { description: "optionalScalarIntegerPrimitiveInput doc" }
    optionalScalarLongInput: { description: "optionalScalarLongInput doc" }
    optionalScalarLongPrimitiveInput: { description: "optionalScalarLongPrimitiveInput doc" }
    optionalScalarStringInput: { description: "optionalScalarStringInput doc" }
  }
}

