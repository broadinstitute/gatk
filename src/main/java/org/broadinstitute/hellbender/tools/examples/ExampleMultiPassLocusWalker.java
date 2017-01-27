package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Example multi-pass locus walker that prints locus-based coverage from supplied reads to the specified output file (stdout if none provided), along with overlapping reference bases/features (if provided)",
        oneLineSummary = "Example tool that prints locus-based coverage with optional contextual data",
        programGroup = ReadProgramGroup.class
)
public class ExampleMultiPassLocusWalker extends MultiPassLocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File output = null;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files", optional = true)
    private List<FeatureInput<VariantContext>> variants;

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = output != null ? new PrintStream(output) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(output, e);
        }
    }

    @Override
    public void apply( String inputLabel, AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        // Get pileup and counts
        ReadPileup pileup = alignmentContext.getBasePileup();
        // print the locus and coverage
        outputStream.printf("In " + inputLabel + ": Current locus %s:%d (coverage=%s)\n", alignmentContext.getContig(),
                alignmentContext.getPosition(), pileup.size());
        // print the reference context if available
        if ( referenceContext.hasBackingDataSource() ) {
            outputStream.println("\tReference base(s): " + new String(referenceContext.getBases()));
        }
        // print the overlapping variants if there are some
        if(featureContext.hasBackingDataSource()) {
            List<VariantContext> vars = featureContext.getValues(variants);
            if(!vars.isEmpty()) {
                outputStream.println("\tOverlapping variant(s):");
                for (VariantContext variant : vars) {
                    outputStream.printf("\t\t%s:%d-%d, Ref:%s, Alt(s):%s\n", variant.getContig(), variant.getStart(),
                            variant.getEnd(), variant.getReference(), variant.getAlternateAlleles());
                }
            }
        }
    }

    @Override
    public void closeTool() {
        if ( outputStream != null )
            outputStream.close();
    }
}
