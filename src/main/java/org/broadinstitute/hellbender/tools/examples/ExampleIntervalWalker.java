package org.broadinstitute.hellbender.tools.examples;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.IntervalProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Example/toy program that shows how to implement the IntervalWalker interface. Prints supplied intervals
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        usage = "Prints intervals supplied via -L to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        usageShort = "Print intervals with optional contextual data",
        programGroup = IntervalProgramGroup.class
)
public class ExampleIntervalWalker extends IntervalWalker {

    @ArgumentCollection
    public OptionalVariantInputArgumentCollection optionalVariants = new OptionalVariantInputArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    public File outputFile = null;

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }
    }

    @Override
    public void apply( SimpleInterval interval, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
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
        for ( SAMRecord read : readsContext ) {
            outputStream.printf("\tOverlapping read at %s:%d-%d\n", read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
        }
        outputStream.println();
    }

    private void printVariants( final FeatureContext featureContext ) {
        for ( VariantContext variant : featureContext.getValues(optionalVariants.variantFiles) ) {
            outputStream.printf("\tOverlapping variant at %s:%d-%d. Ref: %s Alt(s): %s\n",
                    variant.getContig(), variant.getStart(), variant.getEnd(), variant.getReference(), variant.getAlternateAlleles());
        }
        outputStream.println();
    }

    @Override
    public Object onTraversalDone() {
        if ( outputStream != null )
            outputStream.close();

        return null;
    }

}
