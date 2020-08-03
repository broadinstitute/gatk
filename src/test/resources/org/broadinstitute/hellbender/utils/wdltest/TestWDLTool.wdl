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

  }

  output {
    # Workflow Outputs                                  
    File TestWDLToolrequiredScalarFileOutput = TestWDLTool.TestWDLTool_requiredScalarFileOutput
    File TestWDLToolrequiredScalarFileOutputDictionary = TestWDLTool.TestWDLTool_requiredScalarFileOutputDictionary
    File TestWDLToolrequiredScalarFileOutputIndex = TestWDLTool.TestWDLTool_requiredScalarFileOutputIndex
    Array[File] TestWDLToolrequiredListFileOutput = TestWDLTool.TestWDLTool_requiredListFileOutput
    Array[File] TestWDLToolrequiredListFileOutputDictionary = TestWDLTool.TestWDLTool_requiredListFileOutputDictionary
    Array[File] TestWDLToolrequiredListFileOutputIndex = TestWDLTool.TestWDLTool_requiredListFileOutputIndex
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

  }

  command <<<
    ~{gatk} TestWDLTool \
    ~{sep=' ' positionalArgs} \
    --requiredListFileInput ~{sep=' --requiredListFileInput ' requiredListFileInput} \
    --requiredListFileOutput ~{sep=' --requiredListFileOutput ' requiredListFileOutput} \
    --requiredScalarFileInput ~{sep=' --requiredScalarFileInput ' requiredScalarFileInput} \
    --requiredScalarFileOutput ~{sep=' --requiredScalarFileOutput ' requiredScalarFileOutput} \

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
    File TestWDLTool_requiredScalarFileOutput = "${requiredScalarFileOutput}"
    File TestWDLTool_requiredScalarFileOutputDictionary = "${requiredScalarFileOutputDictionary}"
    File TestWDLTool_requiredScalarFileOutputIndex = "${requiredScalarFileOutputIndex}"
    Array[File] TestWDLTool_requiredListFileOutput = "${requiredListFileOutput}"
    Array[File] TestWDLTool_requiredListFileOutputDictionary = "${requiredListFileOutputDictionary}"
    Array[File] TestWDLTool_requiredListFileOutputIndex = "${requiredListFileOutputIndex}"
  }
 }

