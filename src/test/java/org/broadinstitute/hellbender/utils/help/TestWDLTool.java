//package org.broadinstitute.hellbender.utils.help;
//
//import org.broadinstitute.barclay.argparser.Argument;
//import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
//import org.broadinstitute.barclay.argparser.PositionalArguments;
//import org.broadinstitute.barclay.argparser.WorkflowProperties;
//import org.broadinstitute.barclay.argparser.WorkflowInput;
//import org.broadinstitute.barclay.argparser.WorkflowOutput;
//import org.broadinstitute.barclay.help.DocumentedFeature;
//import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
//
//import java.io.File;
//import java.util.List;
//
///**
// * NOTE: this file needs to live in a separate package from the doc tests, otherwise all of the docs tests
// * will pick it up as a command line program, which will change the outputs.
// *
// * CommandLineProgram test tool for testing WDL generation. Contains various combinations of
// * commandline argument and workflow input/outputs with companion resources:
// *
// *  scalar/array
// *  file/non-file
// *  required/optional
// */
//@CommandLineProgramProperties(
//        summary = TestWDLTool.SUMMARY,
//        oneLineSummary = TestWDLTool.ONE_LINE_SUMMARY,
//        programGroup = TestProgramGroup.class)
//@WorkflowProperties(memory ="8G")
//@DocumentedFeature(groupName = TestWDLTool.GROUP_NAME)
//public class TestWDLTool {
//
//    public static final String SUMMARY = "WDL Test Tool";
//    public static final String ONE_LINE_SUMMARY = "WDL Test Tool to test WDL Generation";
//    public static final String GROUP_NAME = "WDL feature group name";
//
//    @PositionalArguments(doc = "Positional args doc")
//    @WorkflowInput(requiredCompanions={"posDictionary", "posIndex"})
//    public List<File> positionalListFileInput;
//
//    // required input Files
//
//    @Argument(fullName = "requiredScalarFileInputNoCompanions",
//            shortName = "requiredScalarFileInputNoCompanions",
//            doc = "requiredScalarFileInputNoCompanions doc",
//            optional = false)
//    @WorkflowInput
//    public File requiredScalarFileInputNoCompanions;
//
//    @Argument(fullName = "requiredScalarFileInputRequiredCompanions",
//            shortName = "requiredScalarFileInputRequiredCompanions",
//            doc = "requiredScalarFileInputRequiredCompanions doc",
//            optional = false)
//    @WorkflowInput(requiredCompanions={"requiredScalarFileInputRequiredCompanionsDictionary", "requiredScalarFileInputRequiredCompanionsIndex"}, localizationOptional = true)
//    public File requiredScalarFileInputRequiredCompanions;
//
//    @Argument(fullName = "requiredScalarFileInputOptionalCompanions",
//            shortName = "requiredScalarFileInputOptionalCompanions",
//            doc = "requiredScalarFileInputOptionalCompanions doc",
//            optional = false)
//    @WorkflowInput(optionalCompanions={"requiredScalarFileInputOptionalCompanionsDictionary", "requiredScalarFileInputOptionalCompanionsIndex"}, localizationOptional = true)
//    public File requiredScalarFileInputOptionalCompanions;
//
//    @Argument(fullName = "requiredListFileInputNoCompanions",
//            shortName = "requiredListFileInputNoCompanions",
//            doc = "requiredListFileInputNoCompanions doc",
//            optional = false)
//    @WorkflowInput
//    public List<File> requiredListFileInputNoCompanions;
//
//    @Argument(fullName = "requiredListFileInputRequiredCompanions",
//            shortName = "requiredListFileInputRequiredCompanions",
//            doc = "requiredListFileInputRequiredCompanions doc",
//            optional = false)
//    @WorkflowInput(requiredCompanions={"requiredListFileInputRequiredCompanionsDictionary", "requiredListFileInputRequiredCompanionsIndex"}, localizationOptional = true)
//    public List<File> requiredListFileInputRequiredCompanions;
//
//    @Argument(fullName = "requiredListFileInputOptionalCompanions",
//            shortName = "requiredListFileInputOptionalCompanions",
//            doc = "requiredListFileInputOptionalCompanions doc",
//            optional = false)
//    @WorkflowInput(optionalCompanions={"requiredListFileInputOptionalCompanionsDictionary", "requiredListFileInputOptionalCompanionsIndex"}, localizationOptional = true)
//    public List<File> requiredListFileInputOptionalCompanions;
//
//    @Argument(fullName = "requiredListFileInputMixedCompanions",
//            shortName = "requiredListFileInputMixedCompanions",
//            doc = "requiredListFileInputMixedCompanions doc",
//            optional = false)
//    @WorkflowInput(
//            requiredCompanions = {"requiredListFileInputMixedCompanionsRequired"},
//            optionalCompanions = {"requiredListFileInputMixedCompanionsOptional"})
//    public List<File> requiredListFileInputMixedCompanions;
//
//    // required output Files
//
//    @Argument(fullName = "requiredScalarFileOutputNoCompanions",
//            shortName = "requiredScalarFileOutputNoCompanions",
//            doc = "requiredScalarFileOutputNoCompanions doc",
//            optional = false)
//    @WorkflowOutput
//    public File requiredScalarFileOutputNoCompanions;
//
//    @Argument(fullName = "requiredScalarFileOutputRequiredCompanions",
//            shortName = "requiredScalarFileOutputRequiredCompanions",
//            doc = "requiredScalarFileOutputRequiredCompanions doc",
//            optional = false)
//    @WorkflowOutput(requiredCompanions={"requiredScalarFileOutputRequiredCompanionsDictionary", "requiredScalarFileOutputRequiredCompanionsIndex"})
//    public File requiredScalarFileOutputRequiredCompanions;
//
//    @Argument(fullName = "requiredScalarFileOutputOptionalCompanions",
//            shortName = "requiredScalarFileOutputOptionalCompanions",
//            doc = "requiredScalarFileOutputOptionalCompanions doc",
//            optional = false)
//    @WorkflowOutput(optionalCompanions={"requiredScalarFileOutputOptionalCompanionsDictionary", "requiredScalarFileOutputOptionalCompanionsIndex"})
//    public File requiredScalarFileOutputOptionalCompanions;
//
//    @Argument(fullName = "requiredListFileOutputNoCompanions",
//            shortName = "requiredListFileOutputNoCompanions",
//            doc = "requiredListFileOutputNoCompanions doc",
//            optional = false)
//    @WorkflowOutput
//    public List<File> requiredListFileOutputNoCompanions;
//
//    @Argument(fullName = "requiredListFileOutputRequiredCompanions",
//            shortName = "requiredListFileOutputRequiredCompanions",
//            doc = "requiredListFileOutputRequiredCompanions doc",
//            optional = false)
//    @WorkflowOutput(requiredCompanions={"requiredListFileOutputRequiredCompanionsDictionary", "requiredListFileOutputRequiredCompanionsIndex"})
//    public List<File> requiredListFileOutputRequiredCompanions;
//
//    @Argument(fullName = "requiredListFileOutputOptionalCompanions",
//            shortName = "requiredListFileOutputOptionalCompanions",
//            doc = "requiredListFileOutputOptionalCompanions doc",
//            optional = false)
//    @WorkflowOutput(optionalCompanions={"requiredListFileOutputOptionalCompanionsDictionary", "requiredListFileOutputOptionalCompanionsIndex"})
//    public List<File> requiredListFileOutputOptionalCompanions;
//
//    @Argument(fullName = "requiredListFileOutputMixedCompanions",
//            shortName = "requiredListFileOutputMixedCompanions",
//            doc = "requiredListFileOutputMixedCompanions doc",
//            optional = false)
//    @WorkflowOutput(
//            requiredCompanions = {"requiredListFileOutputMixedCompanionsRequired"},
//            optionalCompanions = {"requiredListFileOutputMixedCompanionsOptional"})
//    public List<File> requiredListFileOutputMixedCompanions;
//
//    // optional input Files
//
//    @Argument(fullName = "optionalScalarFileInputNoCompanions",
//            shortName = "optionalScalarFileInputNoCompanions",
//            doc = "optionalScalarFileInputNoCompanions doc",
//            optional = true)
//    @WorkflowInput
//    public File optionalScalarFileInputNoCompanions;
//
//    @Argument(fullName = "optionalScalarFileInputOptionalCompanions",
//            shortName = "optionalScalarFileInputOptionalCompanions",
//            doc = "optionalScalarFileInputOptionalCompanions doc",
//            optional = true)
//    @WorkflowInput(optionalCompanions={"optionalScalarFileInputOptionalCompanionsDictionary", "optionalScalarFileInputOptionalCompanionsIndex"})
//    public File optionalScalarFileInputOptionalCompanions;
//
//    @Argument(fullName = "optionalScalarFileInputRequiredCompanions",
//            shortName = "optionalScalarFileInputRequiredCompanions",
//            doc = "optionalScalarFileInputRequiredCompanions doc",
//            optional = true)
//    @WorkflowInput(requiredCompanions={"optionalScalarFileInputRequiredCompanionsDictionary", "optionalScalarFileInputRequiredCompanionsIndex"})
//    public File optionalScalarFileInputRequiredCompanions;
//
//    @Argument(fullName = "optionalListFileInputNoCompanions",
//            shortName = "optionalListFileInputNoCompanions",
//            doc = "optionalListFileInputNoCompanions doc",
//            optional = true)
//    @WorkflowInput
//    public List<File> optionalListFileInputNoCompanions;
//
//    @Argument(fullName = "optionalListFileInputRequiredCompanions",
//            shortName = "optionalListFileInputRequiredCompanions",
//            doc = "optionalListFileInputRequiredCompanions doc",
//            optional = true)
//    @WorkflowInput(requiredCompanions={"optionalListFileInputRequiredCompanionsDictionary", "optionalListFileInputRequiredCompanionsIndex"})
//    public List<File> optionalListFileInputRequiredCompanions;
//
//    @Argument(fullName = "optionalListFileInputOptionalCompanions",
//            shortName = "optionalListFileInputOptionalCompanions",
//            doc = "optionalListFileInputOptionalCompanions doc",
//            optional = true)
//    @WorkflowInput(optionalCompanions={"optionalListFileInputOptionalCompanionsDictionary", "optionalListFileInputOptionalCompanionsIndex"})
//    public List<File> optionalListFileInputOptionalCompanions;
//
//    @Argument(fullName = "optionalListFileInputMixedCompanions",
//            shortName = "optionalListFileInputMixedCompanions",
//            doc = "optionalListFileInputMixedCompanions doc",
//            optional = true)
//    @WorkflowInput(
//            requiredCompanions = {"optionalListFileInputMixedCompanionsRequired"},
//            optionalCompanions = {"optionalListFileInputMixedCompanionsOptional"})
//    public List<File> optionalListFileInputMixedCompanions;
//
//    // optional output Files
//
//    @Argument(fullName = "optionalScalarFileOutputNoCompanions",
//            shortName = "optionalScalarFileOutputNoCompanions",
//            doc = "optionalScalarFileOutputNoCompanions doc",
//            optional = true)
//    @WorkflowOutput
//    public File optionalScalarFileOutputNoCompanions;
//
//    @Argument(fullName = "optionalScalarFileOutputRequiredCompanions",
//            shortName = "optionalScalarFileOutputRequiredCompanions",
//            doc = "optionalScalarFileOutputRequiredCompanions doc",
//            optional = true)
//    @WorkflowOutput(requiredCompanions={"optionalScalarFileOutputRequiredCompanionsDictionary", "optionalScalarFileOutputRequiredCompanionsIndex"})
//    public File optionalScalarFileOutputRequiredCompanions;
//
//    @Argument(fullName = "optionalScalarFileOutputOptionalCompanions",
//            shortName = "optionalScalarFileOutputOptionalCompanions",
//            doc = "optionalScalarFileOutputOptionalCompanions doc",
//            optional = true)
//    @WorkflowOutput(requiredCompanions={"optionalScalarFileOutputOptionalCompanionsDictionary", "optionalScalarFileOutputOptionalCompanionsIndex"})
//    public File optionalScalarFileOutputOptionalCompanions;
//
//    @Argument(fullName = "optionalListFileOutputRequiredCompanions",
//            shortName = "optionalListFileOutputRequiredCompanions",
//            doc = "optionalListFileOutputRequiredCompanions doc",
//            optional = true)
//    @WorkflowOutput(requiredCompanions={"optionalListFileOutputRequiredCompanionsDictionary", "optionalListFileOutputRequiredCompanionsIndex"})
//    public List<File> optionalListFileOutputRequiredCompanions;
//
//    @Argument(fullName = "optionalListFileOutputMixedCompanions",
//            shortName = "optionalListFileOutputMixedCompanions",
//            doc = "optionalListFileOutputMixedCompanions doc",
//            optional = true)
//    @WorkflowOutput(
//            requiredCompanions = {"optionalListFileOutputMixedCompanionsRequired"},
//            optionalCompanions = {"optionalListFileOutputMixedCompanionsOptional"})
//    public List<File> optionalListFileOutputMixedCompanions;
//
//    // non-File types
//
//    @Argument(fullName = "optionalScalarStringInput",
//            shortName = "optionalScalarStringInput",
//            doc = "optionalScalarStringInput doc",
//            optional = true)
//    public String optionalScalarStringInput;
//
//    @Argument(fullName = "optionalListStringInput",
//            shortName = "optionalListStringInput",
//            doc = "optionalListStringInput doc",
//            optional = true)
//    public List<String> optionalListStringInput;
//
//    @Argument(fullName = "optionalScalarIntegerPrimitiveInput",
//            shortName = "optionalScalarIntegerPrimitiveInput",
//            doc = "optionalScalarIntegerPrimitiveInput doc",
//            optional = true)
//    public int optionalScalarIntegerPrimitiveInput;
//
//    @Argument(fullName = "optionalScalarIntegerInput",
//            shortName = "optionalScalarIntegerInput",
//            doc = "optionalScalarIntegerInput doc",
//            optional = true)
//    public Integer optionalScalarIntegerInput;
//
//    @Argument(fullName = "optionalListIntegerInput",
//            shortName = "optionalListIntegerInput",
//            doc = "optionalListIntegerInput doc",
//            optional = true)
//    public List<Integer> optionalListIntegerInput;
//
//    @Argument(fullName = "optionalScalarLongPrimitiveInput",
//            shortName = "optionalScalarLongPrimitiveInput",
//            doc = "optionalScalarLongPrimitiveInput doc",
//            optional = true)
//    public long optionalScalarLongPrimitiveInput;
//
//    @Argument(fullName = "optionalScalarLongInput",
//            shortName = "optionalScalarLongInput",
//            doc = "optionalScalarLongInput doc",
//            optional = true)
//    public Long optionalScalarLongInput;
//
//    @Argument(fullName = "optionalListLongInput",
//            shortName = "optionalListLongInput",
//            doc = "optionalListLongInput doc",
//            optional = true)
//    public List<Long> optionalListLongInput;
//
//    @Argument(fullName = "optionalScalarFloatPrimitiveInput",
//            shortName = "optionalScalarFloatPrimitiveInput",
//            doc = "optionalScalarFloatPrimitiveInput doc",
//            optional = true)
//    public float optionalScalarFloatPrimitiveInput;
//
//    @Argument(fullName = "optionalScalarFloatInput",
//            shortName = "optionalScalarFloatInput",
//            doc = "optionalScalarFloatInput doc",
//            optional = true)
//    public Float optionalScalarFloatInput;
//
//    @Argument(fullName = "optionalListFloatInput",
//            shortName = "optionalListFloatInput",
//            doc = "optionalListFloatInput doc",
//            optional = true)
//    public List<Float> optionalListFloatInput;
//
//    @Argument(fullName = "optionalScalarDoublePrimitiveInput",
//            shortName = "optionalScalarDoublePrimitiveInput",
//            doc = "optionalScalarDoublePrimitiveInput doc",
//            optional = true)
//    public double optionalScalarDoublePrimitiveInput;
//
//    @Argument(fullName = "optionalScalarDoubleInput",
//            shortName = "optionalScalarDoubleInput",
//            doc = "optionalScalarDoubleInput doc",
//            optional = true)
//    public Double optionalScalarDoubleInput;
//
//    @Argument(fullName = "optionalListDoubleInput",
//            shortName = "optionalListDoubleInput",
//            doc = "optionalListDoubleInput doc",
//            optional = true)
//    public List<Double> optionalListDoubleInput;
//
//}
