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

    @Argument(fullName = "architecture", shortName = "a", doc = "Neural Net architecture and weights hd5 file", optional = false)
    private String architecture = null;

    @Argument(fullName = "batchSize", shortName = "bs", doc = "Batch size", optional = true)
    private String batchSize = "32";

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

        System.out.println("output vcf is:"+outputFile.getAbsolutePath());
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        //outputStream.println("Current variant: " + variant);
    }

    @Override
    public Object onTraversalSuccess() {
        outputStream.close();
        outputStream = null;

        final Resource pythonScriptResource = new Resource("ApplyNeuralNetModel.py", NeuralNetExecutor.class);
        System.out.println("Run neural Net script");
        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                Arrays.asList(
                        "--architecture",
                        architecture,
                        "--input_vcf",
                        drivingVariantFile,
                        "--output_vcf",
                        outputFile.getAbsolutePath(),
                        "--reference_fasta",
                        referenceArguments.getReferenceFileName(),
                        "--batch_size",
                        batchSize
                )

        );
        System.out.println("Ran neural Net script");
        return pythonReturnCode;
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
