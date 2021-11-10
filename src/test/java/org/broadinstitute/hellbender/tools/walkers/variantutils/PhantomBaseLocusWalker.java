package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.PrintStream;
import java.util.List;
import java.util.stream.Collectors;


@DocumentedFeature
@CommandLineProgramProperties(
        oneLineSummary = "Walk through reads at site and look at phantom bases in hmers.",
        summary = "Walk bam and look at phantom bases in hmers.",
        programGroup = ReferenceProgramGroup.class
)

public class PhantomBaseLocusWalker extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private GATKPath outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    @Argument(shortName = "F", fullName = "feature_file", doc = "Feature file (eg., VCF or BED file)")
    public FeatureInput<TableFeature> featuresFile = null;

    private PrintStream outputStream = null;



    @Override
    public void onTraversalStart() {
        outputStream = outputFile != null ? new PrintStream(outputFile.getOutputStream()) : System.out;
//        outputStream.println("HEADER GT GQ MaxHmerLength Call");
    }




    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        final String test = "Hello World";

        if ( ! featureContext.getValues(featuresFile).isEmpty() ) {
            List<GATKRead> highMapQReads = alignmentContext.getBasePileup().getReads().stream().filter(read -> read.getMappingQuality() >= 30).collect(Collectors.toList());


        }

        // alignmentContext.getBasePileup().getElementStream().filter(pe -> pe.getRead().getMappingQuality() > 50).collect(Collectors.toList())

    }


    @Override
    public void closeTool() {


        if ( outputStream != null ) {
            outputStream.close();
        }
    }




}
