version 1.0

# Run TestWDLTool (WDL auto generated from: GATK Version 1.1-111)
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
#    posDictionary                                      Companion resource for: Positional args doc                 
#    posIndex                                           Companion resource for: Positional args doc                 
#
#  Required Tool Arguments
#    requiredListFileInput                              requiredListFileInput doc                                   
#    requiredListFileInputDictionary                    Companion resource for: requiredListFileInput doc           
#    requiredListFileInputIndex                         Companion resource for: requiredListFileInput doc           
#    requiredListFileOutput                             requiredListFileOutput doc                                  
#    requiredListFileOutputDictionary                   Companion resource for: requiredListFileOutput doc          
#    requiredListFileOutputIndex                        Companion resource for: requiredListFileOutput doc          
#    requiredScalarFileInput                            requiredScalarFileInput doc                                 
#    requiredScalarFileInputDictionary                  Companion resource for: requiredScalarFileInput doc         
#    requiredScalarFileInputIndex                       Companion resource for: requiredScalarFileInput doc         
#    requiredScalarFileOutput                           requiredScalarFileOutput doc                                
#    requiredScalarFileOutputDictionary                 Companion resource for: requiredScalarFileOutput doc        
#    requiredScalarFileOutputIndex                      Companion resource for: requiredScalarFileOutput doc        
#
#  Optional Tool Arguments
#    optionaldListFileOutput                            optionalListFileOutput doc                                  
#    optionalListFileOutputDictionary                   Companion resource for: optionalListFileOutput doc          
#    requiredListFileOutputIndex                        Companion resource for: optionalListFileOutput doc          
#    optionaldScalarFileOutput                          optionalScalarFileOutput doc                                
#    optionalScalarFileOutputDictionary                 Companion resource for: optionalScalarFileOutput doc        
#    requiredScalarFileOutputIndex                      Companion resource for: optionalScalarFileOutput doc        
#    optionalListDoubleInput                            optionalListDoubleInput doc                                 
#    optionalListFileInput                              optionalListFileInput doc                                   
#    optionalListFileInputDictionary                    Companion resource for: optionalListFileInput doc           
#    requiredListFileInputIndex                         Companion resource for: optionalListFileInput doc           
#    optionalListFloatInput                             optionalListFloatInput doc                                  
#    optionalListIntegerInput                           optionalListIntegerInput doc                                
#    optionalListLongInput                              optionalListLongInput doc                                   
#    optionalListStringInput                            optionalListStringInput doc                                 
#    optionalScalarDoubleInput                          optionalScalarDoubleInput doc                               
#    optionalScalarDoublePrimitiveInput                 optionalScalarDoublePrimitiveInput doc                      
#    optionalScalarFileInput                            optionalScalarFileInput doc                                 
#    optionalScalarFileInputDictionary                  Companion resource for: optionalScalarFileInput doc         
#    requiredScalarFileInputIndex                       Companion resource for: optionalScalarFileInput doc         
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
    Array[File] requiredListFileInput
    Array[File] requiredListFileInputDictionary
    Array[File] requiredListFileInputIndex
    Array[String] requiredListFileOutput
    Array[String] requiredListFileOutputDictionary
    Array[String] requiredListFileOutputIndex
    File requiredScalarFileInput
    File requiredScalarFileInputDictionary
    File requiredScalarFileInputIndex
    String requiredScalarFileOutput
    String requiredScalarFileOutputDictionary
    String requiredScalarFileOutputIndex

    # Optional Tool Arguments
    Array[String]? optionaldListFileOutput
    Array[String]? optionalListFileOutputDictionary
    Array[String]? requiredListFileOutputIndex
    String? optionaldScalarFileOutput
    String? optionalScalarFileOutputDictionary
    String? requiredScalarFileOutputIndex
    Array[Float]? optionalListDoubleInput
    Array[File]? optionalListFileInput
    Array[File]? optionalListFileInputDictionary
    Array[File]? requiredListFileInputIndex
    Array[Float]? optionalListFloatInput
    Array[Int]? optionalListIntegerInput
    Array[Int]? optionalListLongInput
    Array[String]? optionalListStringInput
    Float? optionalScalarDoubleInput
    Float? optionalScalarDoublePrimitiveInput
    File? optionalScalarFileInput
    File? optionalScalarFileInputDictionary
    File? requiredScalarFileInputIndex
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
        requiredListFileInput                              = requiredListFileInput,
        requiredListFileInputDictionary                    = requiredListFileInputDictionary,
        requiredListFileInputIndex                         = requiredListFileInputIndex,
        requiredListFileOutput                             = requiredListFileOutput,
        requiredListFileOutputDictionary                   = requiredListFileOutputDictionary,
        requiredListFileOutputIndex                        = requiredListFileOutputIndex,
        requiredScalarFileInput                            = requiredScalarFileInput,
        requiredScalarFileInputDictionary                  = requiredScalarFileInputDictionary,
        requiredScalarFileInputIndex                       = requiredScalarFileInputIndex,
        requiredScalarFileOutput                           = requiredScalarFileOutput,
        requiredScalarFileOutputDictionary                 = requiredScalarFileOutputDictionary,
        requiredScalarFileOutputIndex                      = requiredScalarFileOutputIndex,

        # Optional Tool Arguments
        optionaldListFileOutput                            = optionaldListFileOutput,
        optionalListFileOutputDictionary                   = optionalListFileOutputDictionary,
        requiredListFileOutputIndex                        = requiredListFileOutputIndex,
        optionaldScalarFileOutput                          = optionaldScalarFileOutput,
        optionalScalarFileOutputDictionary                 = optionalScalarFileOutputDictionary,
        requiredScalarFileOutputIndex                      = requiredScalarFileOutputIndex,
        optionalListDoubleInput                            = optionalListDoubleInput,
        optionalListFileInput                              = optionalListFileInput,
        optionalListFileInputDictionary                    = optionalListFileInputDictionary,
        requiredListFileInputIndex                         = requiredListFileInputIndex,
        optionalListFloatInput                             = optionalListFloatInput,
        optionalListIntegerInput                           = optionalListIntegerInput,
        optionalListLongInput                              = optionalListLongInput,
        optionalListStringInput                            = optionalListStringInput,
        optionalScalarDoubleInput                          = optionalScalarDoubleInput,
        optionalScalarDoublePrimitiveInput                 = optionalScalarDoublePrimitiveInput,
        optionalScalarFileInput                            = optionalScalarFileInput,
        optionalScalarFileInputDictionary                  = optionalScalarFileInputDictionary,
        requiredScalarFileInputIndex                       = requiredScalarFileInputIndex,
        optionalScalarFloatInput                           = optionalScalarFloatInput,
        optionalScalarFloatPrimitiveInput                  = optionalScalarFloatPrimitiveInput,
        optionalScalarIntegerInput                         = optionalScalarIntegerInput,
        optionalScalarIntegerPrimitiveInput                = optionalScalarIntegerPrimitiveInput,
        optionalScalarLongInput                            = optionalScalarLongInput,
        optionalScalarLongPrimitiveInput                   = optionalScalarLongPrimitiveInput,
        optionalScalarStringInput                          = optionalScalarStringInput,

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
    Array[File] requiredListFileInput
    Array[File] requiredListFileInputDictionary
    Array[File] requiredListFileInputIndex
    Array[String] requiredListFileOutput
    Array[String] requiredListFileOutputDictionary
    Array[String] requiredListFileOutputIndex
    File requiredScalarFileInput
    File requiredScalarFileInputDictionary
    File requiredScalarFileInputIndex
    String requiredScalarFileOutput
    String requiredScalarFileOutputDictionary
    String requiredScalarFileOutputIndex
    Array[String]? optionaldListFileOutput
    Array[String]? optionalListFileOutputDictionary
    Array[String]? requiredListFileOutputIndex
    String? optionaldScalarFileOutput
    String? optionalScalarFileOutputDictionary
    String? requiredScalarFileOutputIndex
    Array[Float]? optionalListDoubleInput
    Array[File]? optionalListFileInput
    Array[File]? optionalListFileInputDictionary
    Array[File]? requiredListFileInputIndex
    Array[Float]? optionalListFloatInput
    Array[Int]? optionalListIntegerInput
    Array[Int]? optionalListLongInput
    Array[String]? optionalListStringInput
    Float? optionalScalarDoubleInput
    Float? optionalScalarDoublePrimitiveInput
    File? optionalScalarFileInput
    File? optionalScalarFileInputDictionary
    File? requiredScalarFileInputIndex
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
    --requiredListFileInput ~{sep=' --requiredListFileInput ' requiredListFileInput} \
    --requiredListFileOutput ~{sep=' --requiredListFileOutput ' requiredListFileOutput} \
    --requiredScalarFileInput ~{sep=' --requiredScalarFileInput ' requiredScalarFileInput} \
    --requiredScalarFileOutput ~{sep=' --requiredScalarFileOutput ' requiredScalarFileOutput} \
    ~{true='--optionaldListFileOutput ' false='' defined(optionaldListFileOutput)}~{sep=' --optionaldListFileOutput ' optionaldListFileOutput} \
    ~{true='--optionaldScalarFileOutput ' false='' defined(optionaldScalarFileOutput)}~{sep=' --optionaldScalarFileOutput ' optionaldScalarFileOutput} \
    ~{true='--optionalListDoubleInput ' false='' defined(optionalListDoubleInput)}~{sep=' --optionalListDoubleInput ' optionalListDoubleInput} \
    ~{true='--optionalListFileInput ' false='' defined(optionalListFileInput)}~{sep=' --optionalListFileInput ' optionalListFileInput} \
    ~{true='--optionalListFloatInput ' false='' defined(optionalListFloatInput)}~{sep=' --optionalListFloatInput ' optionalListFloatInput} \
    ~{true='--optionalListIntegerInput ' false='' defined(optionalListIntegerInput)}~{sep=' --optionalListIntegerInput ' optionalListIntegerInput} \
    ~{true='--optionalListLongInput ' false='' defined(optionalListLongInput)}~{sep=' --optionalListLongInput ' optionalListLongInput} \
    ~{true='--optionalListStringInput ' false='' defined(optionalListStringInput)}~{sep=' --optionalListStringInput ' optionalListStringInput} \
    ~{true='--optionalScalarDoubleInput ' false='' defined(optionalScalarDoubleInput)}~{sep=' --optionalScalarDoubleInput ' optionalScalarDoubleInput} \
    ~{true='--optionalScalarDoublePrimitiveInput ' false='' defined(optionalScalarDoublePrimitiveInput)}~{sep=' --optionalScalarDoublePrimitiveInput ' optionalScalarDoublePrimitiveInput} \
    ~{true='--optionalScalarFileInput ' false='' defined(optionalScalarFileInput)}~{sep=' --optionalScalarFileInput ' optionalScalarFileInput} \
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

 }

