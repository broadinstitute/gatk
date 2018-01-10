package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Example/toy program that shows how to implement the AssemblyRegionWalker interface.
 *
 * Prints out the bounds of each assembly region with and without padding, as well as the number of reads in each region.
 * Also prints overlapping reference bases / variants for each region, if provided.
 */
@CommandLineProgramProperties(
        summary = "Example AssemblyRegionWalker that prints out the bounds of each assembly region with and without padding, as well as the number of reads in each region",
        oneLineSummary = "Example AssemblyRegionWalker",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleAssemblyRegionWalker extends AssemblyRegionWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

    @Argument(fullName="knownVariants", shortName="knownVariants", doc="Known set of variants", optional=true)
    private FeatureInput<VariantContext> knownVariants;

    private PrintStream outputStream = null;

    @Override
    protected int defaultMinAssemblyRegionSize() { return 50; }

    @Override
    protected int defaultMaxAssemblyRegionSize() { return 300; }

    @Override
    protected int defaultAssemblyRegionPadding() { return 100; }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() { return 50; }

    @Override
    protected double defaultActiveProbThreshold() { return 0.002; }

    @Override
    protected int defaultMaxProbPropagationDistance() { return 50; }

    @Override
    protected boolean includeReadsWithDeletionsInIsActivePileups() { return true; }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() {
        // This example AssemblyRegionEvaluator considers all loci to be active
        // (ie., it assigns an isActive probability of 1.0 to all loci):
        return (locusPileup, referenceContext, featureContext) -> new ActivityProfileState(new SimpleInterval(locusPileup), 1.0);
    }

    @Override
    public void onTraversalStart() {
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public void apply( AssemblyRegion region, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputStream.printf("%s assembly region at %s (%s with padding), containing %d reads.\n\n",
                region.isActive() ? "ACTIVE" : "INACTIVE", region.getSpan(), region.getExtendedSpan(), region.getReads().size());

        printReferenceBases(referenceContext);

        if ( featureContext.hasBackingDataSource() ) {
            printOverlappingVariants(featureContext);
        }
    }

    private void printReferenceBases( final ReferenceContext refContext ) {
        outputStream.printf("\tOverlapping reference bases: %s\n\n", new String(refContext.getBases()));
    }

    private void printOverlappingVariants( final FeatureContext featureContext ) {
        for ( final VariantContext variant : featureContext.getValues(knownVariants) ) {
            outputStream.printf("\tOverlapping variant at %s:%d-%d. Ref: %s Alt(s): %s\n\n",
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
