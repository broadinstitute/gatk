package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

/**
 * Example/toy program that shows how to implement the VariantWalker interface. Prints supplied variants
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints variants supplied to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Example tool that prints variants with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class BlahVariantWalker extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    private PrintStream outputStream = null;

    private Integer lastSeenPosition = null;

    private final int GQ_CUTOFF = 60;

    @Override
    public void onTraversalStart() {
        // look for duplicate samples and positions covered
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final String sampleName = variant.getGenotype(0).getSampleName();
        outputStream.println("Current variant: " + variant);
        // This is a wrapper to loop thru createTSV function -- and split out the VET and PET tables

        // create VET output
        if (variant.isVariant()) {
            try {
                final List<String> TSVtoCreate = BlahVetCreation.createTSV(variant);
            } catch (final IOException e) {
                throw new IllegalArgumentException("Current variant is missing required fields", e);
            }
        }
        // create PET output
        try {
            if (variant.getGenotype(0).getGQ() < GQ_CUTOFF) {
                final List<String> TSVtoCreate = BlahPetCreation.createTSV(variant);
            }
        } catch (final Exception e) {
            throw new IllegalArgumentException("GQ NOT GOOD", e);
        }

        // create "missing" variants

        if (!(lastSeenPosition == null) && !(lastSeenPosition + 1 == variant.getStart())){
            // actually make sure this is a position we call over -- may want to use interval lists (ask David)
            BlahPetCreation.createMissingTSV(lastSeenPosition + 1, variant.getEnd(), sampleName);
        }

        lastSeenPosition = variant.getEnd();
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
