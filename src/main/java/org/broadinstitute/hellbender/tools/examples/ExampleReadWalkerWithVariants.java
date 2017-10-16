package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

/**
 * Example/toy program that prints reads from the provided file or files along with overlapping variants
 * (if a source of variants is provided). Intended to show how to implement the ReadWalker interface.
 */
@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) along with overlapping variants (if a source of variants is provided) to the specified output file (or STDOUT if none specified)",
        oneLineSummary = "Print reads with overlapping variants",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleReadWalkerWithVariants extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files", optional = true)
    private List<FeatureInput<VariantContext>> variants;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile;

    @Argument(fullName = "groupVariantsBySource", shortName = "groupVariantsBySource", doc = "If true, group overlapping variants by their source when outputting them", optional = true)
    private boolean groupVariantsBySource = false;

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputStream.printf("Read at %s:%d-%d:\n%s\n", read.getContig(), read.getStart(), read.getEnd(), read.getBasesString());

        if ( groupVariantsBySource ) {
            // We can keep the variants from each source separate by passing in the FeatureInputs
            // individually to featureContext.getValues()
            for ( FeatureInput<VariantContext> featureSource : variants ) {
                outputStream.println("From source " + featureSource.getName());

                for ( VariantContext variant : featureContext.getValues(featureSource) ) {
                    outputStream.printf("\t");
                    printOverlappingVariant(variant);
                }
            }
        }
        else {
            // Passing in all FeatureInputs at once to featureContext.getValues() lets us get
            // all overlapping variants without regard to source
            for ( VariantContext variant : featureContext.getValues(variants) ) {
                printOverlappingVariant(variant);
            }
        }

        outputStream.println();
    }

    private void printOverlappingVariant( final VariantContext variant ) {
        outputStream.printf("Overlapping variant at %s:%d-%d: Ref: %s Alt(s): %s\n",
                variant.getContig(), variant.getStart(), variant.getEnd(),
                variant.getReference(), variant.getAlternateAlleles());
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}