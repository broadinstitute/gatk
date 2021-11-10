package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.PrintStream;
import java.util.List;
import java.util.stream.Collectors;


@DocumentedFeature
@CommandLineProgramProperties(
        oneLineSummary = "Walk through loci and count coverage.",
        summary = "Walk bam and look at coverage by position.",
        programGroup = ReferenceProgramGroup.class
)

public class CoverageWalker extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private GATKPath outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    @Argument(shortName = "F", fullName = "feature_file", doc = "Feature file (eg., VCF or BED file)", optional=true)
    public FeatureInput<TableFeature> featuresFile = null;

//    @Argument(shortName = "eB", fullName = "ExtraBam", doc = "Other bam to compare to")
//    public FeatureInput<AlignmentContext> extraBam;

    private PrintStream outputStream = null;



    @Override
    public void onTraversalStart() {
        outputStream = outputFile != null ? new PrintStream(outputFile.getOutputStream()) : System.out;
        outputStream.println("Rel_Position" + '\t' + "Coverage");
    }




    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        List<GATKRead> highMapQReadsBam = alignmentContext.getBasePileup().getReads().stream().filter(read -> read.getMappingQuality() >= 10).collect(Collectors.toList());

        final int relPosition = alignmentContext.getStart();
        final int coverage = highMapQReadsBam.size();

        outputStream.println(relPosition + "\t" + coverage);

    }


    @Override
    public void closeTool() {

        if ( outputStream != null ) {
            outputStream.close();
        }
    }




}
