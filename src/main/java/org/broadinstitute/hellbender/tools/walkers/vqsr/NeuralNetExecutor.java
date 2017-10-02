package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * VariantWalker which applies VQSR Neural Net Post Traversal.
 */
@CommandLineProgramProperties(
        summary = "Apply a pre-trained neural net model to a directory of read tensors.",
        oneLineSummary = "Execute Neural Network variant filtration in python",
        programGroup = ExampleProgramGroup.class
)

public final class NeuralNetExecutor extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    // Start the Python executor. This does not actually start the Python process, but fails if python can't be located
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        outputStream.println("Current variant: " + variant);

    }

    @Override
    public Object onTraversalSuccess() {
        // Close the stream, and execute a python script that post-processes the output from the java processing
        // step by just copying it to a separate output file:
        //
        // # Test script that accepts two input arguments that are file paths.
        // # Copies the contents of th first file to the second file:
        // import sys
        // with open(sys.argv[1]) as fin:
        //    lines = fin.readlines()
        //    with open(sys.argv[2], "w") as fout:
        //        fout.writelines(lines)
        //
        outputStream.close();
        outputStream = null;

        final Resource pythonScriptResource = new Resource("ApplyNeuralNetModel.py", NeuralNetExecutor.class);
        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                Arrays.asList(
                        "--architecture",
                        "/dsde/working/sam/palantir_cnn/Analysis/vqsr_cnn/weights/m__base_quality_mode_phot__channels_last_False__id_g94982_no_qual_train2__window_size_128__read_limit_128__random_seed_12878__tensor_map_2d_mapping_quality__mode_ref_read_anno.hd5",
                        "--tensors",
                        "/Users/sam/vqsr_data/tensors/not_indel_subset/"
                )
        );

        return pythonReturnCode;
    }


    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
