package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.PositionalArguments;
import org.broadinstitute.barclay.argparser.RuntimeProperties;
import org.broadinstitute.barclay.argparser.WorkflowResource;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;

import java.io.File;
import java.util.List;

/**
 * NOTE: this file needs to live in a separate package from the doc tests, otherwise all of the docs tests
 * will pick it up as a command line program, which will change the outputs.
 *
 * CommandLineProgram test tool for testing WDL generation. Contains various combinations of
 * commandline argument and workflow input/outputs with companion resources:
 *
 *  scalar/array
 *  file/non-file
 *  required/optional
 */
@CommandLineProgramProperties(
        summary = TestWDLTool.SUMMARY,
        oneLineSummary = TestWDLTool.ONE_LINE_SUMMARY,
        programGroup = TestProgramGroup.class)
@RuntimeProperties(memory ="8G")
@DocumentedFeature(groupName = TestWDLTool.GROUP_NAME)
public class TestWDLTool {

    public static final String SUMMARY = "WDL Test Tool";
    public static final String ONE_LINE_SUMMARY = "WDL Test Tool to test WDL Generation";
    public static final String GROUP_NAME = "WDL feature group name";

    @PositionalArguments(doc = "Positional args doc")
    @WorkflowResource(input=true, output=false, companionResources={"posDictionary", "posIndex"})
    public List<File> positionalListFileInput;

    // required Files
    @Argument(fullName = "requiredScalarFileInput",
            shortName = "requiredScalarFileInput",
            doc = "requiredScalarFileInput doc",
            optional = false)
    @WorkflowResource(input=true, output=false, companionResources={"requiredScalarFileInputDictionary", "requiredScalarFileInputIndex"})
    public File requiredScalarFileInput;

    @Argument(fullName = "requiredListFileInput",
            shortName = "requiredListFileInput",
            doc = "requiredListFileInput doc",
            optional = false)
    @WorkflowResource(input=true, output=false, companionResources={"requiredListFileInputDictionary", "requiredListFileInputIndex"})
    public List<File> requiredListFileInput;

    @Argument(fullName = "requiredScalarFileOutput",
            shortName = "requiredScalarFileOutput",
            doc = "requiredScalarFileOutput doc",
            optional = false)
    @WorkflowResource(input=false, output=true, companionResources={"requiredScalarFileOutputDictionary", "requiredScalarFileOutputIndex"})
    public File requiredScalarFileOutput;

    @Argument(fullName = "requiredListFileOutput",
            shortName = "requiredListFileOutput",
            doc = "requiredListFileOutput doc",
            optional = false)
    @WorkflowResource(input=false, output=true, companionResources={"requiredListFileOutputDictionary", "requiredListFileOutputIndex"})
    public List<File> requiredListFileOutput;


    // optional Files
    @Argument(fullName = "optionalScalarFileInput",
            shortName = "optionalScalarFileInput",
            doc = "optionalScalarFileInput doc",
            optional = true)
    @WorkflowResource(input=true, output=false, companionResources={"optionalScalarFileInputDictionary", "requiredScalarFileInputIndex"})
    public File optionalScalarFileInput;

    @Argument(fullName = "optionalListFileInput",
            shortName = "optionalListFileInput",
            doc = "optionalListFileInput doc",
            optional = true)
    @WorkflowResource(input=true, output=false, companionResources={"optionalListFileInputDictionary", "requiredListFileInputIndex"})
    public List<File> optionalListFileInput;

    @Argument(fullName = "optionaldScalarFileOutput",
            shortName = "optionalScalarFileOutput",
            doc = "optionalScalarFileOutput doc",
            optional = true)
    @WorkflowResource(input=false, output=true, companionResources={"optionalScalarFileOutputDictionary", "requiredScalarFileOutputIndex"})
    public File optionalScalarFileOutput;

    @Argument(fullName = "optionaldListFileOutput",
            shortName = "optionalListFileOutput",
            doc = "optionalListFileOutput doc",
            optional = true)
    @WorkflowResource(input=false, output=true, companionResources={"optionalListFileOutputDictionary", "requiredListFileOutputIndex"})
    public List<File> optionalListFileOutput;

    // non-File types

    @Argument(fullName = "optionalScalarStringInput",
            shortName = "optionalScalarStringInput",
            doc = "optionalScalarStringInput doc",
            optional = true)
    public String optionalScalarStringInput;

    @Argument(fullName = "optionalListStringInput",
            shortName = "optionalListStringInput",
            doc = "optionalListStringInput doc",
            optional = true)
    public List<String> optionalListStringInput;

    @Argument(fullName = "optionalScalarIntegerPrimitiveInput",
            shortName = "optionalScalarIntegerPrimitiveInput",
            doc = "optionalScalarIntegerPrimitiveInput doc",
            optional = true)
    public int optionalScalarIntegerPrimitiveInput;

    @Argument(fullName = "optionalScalarIntegerInput",
            shortName = "optionalScalarIntegerInput",
            doc = "optionalScalarIntegerInput doc",
            optional = true)
    public Integer optionalScalarIntegerInput;

    @Argument(fullName = "optionalListIntegerInput",
            shortName = "optionalListIntegerInput",
            doc = "optionalListIntegerInput doc",
            optional = true)
    public List<Integer> optionalListIntegerInput;

    @Argument(fullName = "optionalScalarLongPrimitiveInput",
            shortName = "optionalScalarLongPrimitiveInput",
            doc = "optionalScalarLongPrimitiveInput doc",
            optional = true)
    public long optionalScalarLongPrimitiveInput;

    @Argument(fullName = "optionalScalarLongInput",
            shortName = "optionalScalarLongInput",
            doc = "optionalScalarLongInput doc",
            optional = true)
    public Long optionalScalarLongInput;

    @Argument(fullName = "optionalListLongInput",
            shortName = "optionalListLongInput",
            doc = "optionalListLongInput doc",
            optional = true)
    public List<Long> optionalListLongInput;

    @Argument(fullName = "optionalScalarFloatPrimitiveInput",
            shortName = "optionalScalarFloatPrimitiveInput",
            doc = "optionalScalarFloatPrimitiveInput doc",
            optional = true)
    public float optionalScalarFloatPrimitiveInput;

    @Argument(fullName = "optionalScalarFloatInput",
            shortName = "optionalScalarFloatInput",
            doc = "optionalScalarFloatInput doc",
            optional = true)
    public Float optionalScalarFloatInput;

    @Argument(fullName = "optionalListFloatInput",
            shortName = "optionalListFloatInput",
            doc = "optionalListFloatInput doc",
            optional = true)
    public List<Float> optionalListFloatInput;

    @Argument(fullName = "optionalScalarDoublePrimitiveInput",
            shortName = "optionalScalarDoublePrimitiveInput",
            doc = "optionalScalarDoublePrimitiveInput doc",
            optional = true)
    public double optionalScalarDoublePrimitiveInput;

    @Argument(fullName = "optionalScalarDoubleInput",
            shortName = "optionalScalarDoubleInput",
            doc = "optionalScalarDoubleInput doc",
            optional = true)
    public Double optionalScalarDoubleInput;

    @Argument(fullName = "optionalListDoubleInput",
            shortName = "optionalListDoubleInput",
            doc = "optionalListDoubleInput doc",
            optional = true)
    public List<Double> optionalListDoubleInput;

}
