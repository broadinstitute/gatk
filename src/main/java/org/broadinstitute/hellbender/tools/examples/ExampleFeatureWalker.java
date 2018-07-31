package org.broadinstitute.hellbender.tools.examples;

import htsjdk.tribble.Feature;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Example/toy program that shows how to implement the FeatureWalker interface. Prints supplied features
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints features supplied to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Example tool that prints features with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleFeatureWalker extends FeatureWalker<BEDFeature> {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    private PrintStream outputStream = null;

    @Argument(shortName = "F", fullName = "feature_file", doc = "Feature file (eg., VCF or BED file)")
    public File featuresFile;

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
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(BEDFeature.class);
    }

    @Override
    public File getDrivingFeatureFile() {
        return featuresFile;
    }

    @Override
    public void apply(final BEDFeature feature, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        outputStream.println("Current feature: " + toBedFeatureString(feature));

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

    private static String toBedFeatureString(Feature feature) {
        //Note: BEDFeatures have no toString methods and render as useless pointers
        return feature.getClass().getCanonicalName()+ ":" + feature.getContig() + ":" + feature.getStart() + "-" + feature.getEnd();
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
        for ( final VariantContext variant : featureContext.getValues(auxiliaryVariants) ) {
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
