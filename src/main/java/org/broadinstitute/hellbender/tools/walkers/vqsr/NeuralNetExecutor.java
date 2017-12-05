package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.avro.generic.GenericData;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

    @Argument(fullName = "bamFile", shortName = "bam", doc = "BAM or SAM file", optional = true)
    private String bamFile = null;

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

        logger.info("Output vcf:"+outputFile.getAbsolutePath()+"\nIntervals:"+intervalArgumentCollection.getIntervalStrings().get(0));
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
        List<String> args = new ArrayList<String>();
        if (bamFile != null) {
            args.add("2d");
        } else {
            args.add("1d");
        }
        args.add("--architecture");
        args.add(architecture);
        args.add("--input_vcf");
        args.add(drivingVariantFile);
        args.add("--output_vcf");
        args.add(outputFile.getAbsolutePath());
        args.add("--reference_fasta");
        args.add(referenceArguments.getReferenceFileName());
        args.add("--batch_size");
        args.add(batchSize);
        args.add("--interval_list");
        args.add(intervalArgumentCollection.getIntervalStrings().get(0));

        if (bamFile != null) {
            args.add("--bam_file");
            args.add(bamFile);
        }

        logger.info("Run neural net script");
        final boolean pythonReturnCode = pythonExecutor.executeScript(pythonScriptResource, null, args);
        logger.info("Ran neural net script");
        return pythonReturnCode;
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
