package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Example/toy program that shows how to implement the IntervalWalker interface. Prints supplied intervals
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints intervals supplied via -L to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Print intervals with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleIntervalWalker extends IntervalWalker {

    @ArgumentCollection
    private OptionalVariantInputArgumentCollection optionalVariants = new OptionalVariantInputArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

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
    public void apply( final SimpleInterval interval, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        outputStream.println("Current interval: " + interval);

        if ( referenceContext.hasBackingDataSource() ) {
            printReferenceBases(referenceContext);
        }

        if ( readsContext.hasBackingDataSource() ) {
            printReads(readsContext);
        }

        if ( featureContext.hasBackingDataSource() ) {
            printVariants(featureContext);
        }
    }

    private void printReferenceBases( final ReferenceContext refContext ) {
        outputStream.printf("\tOverlapping reference bases: %s\n\n", new String(refContext.getBases()));
    }

    private void printReads( final ReadsContext readsContext ) {
        for ( final GATKRead read : readsContext ) {
            outputStream.printf("\tOverlapping read at %s:%d-%d\n", read.getContig(), read.getStart(), read.getEnd());
        }
        outputStream.println();
    }

    private void printVariants( final FeatureContext featureContext ) {
        for ( final VariantContext variant : featureContext.getValues(optionalVariants.variantFiles) ) {
            outputStream.printf("\tOverlapping variant at %s:%d-%d. Ref: %s Alt(s): %s\n",
                    variant.getContig(), variant.getStart(), variant.getEnd(), variant.getReference(), variant.getAlternateAlleles());
        }
        outputStream.println();
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }

}
