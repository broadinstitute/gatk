package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.spark.LocusWalkerContext;
import org.broadinstitute.hellbender.engine.spark.LocusWalkerSpark;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.PrintStream;
import java.util.List;

/**
 * Example/toy program that shows how to implement the LocusWalker interface. Prints locus-based coverage from supplied
 * reads, and reference bases/overlapping variants if provided
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints locus-based coverage from supplied read to the specified output file (stdout if none provided), along with overlapping reference bases/features (if provided)",
        oneLineSummary = "Example tool that prints locus-based coverage with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleLocusWalkerSpark extends LocusWalkerSpark {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files", optional = true)
    private List<FeatureInput<VariantContext>> variants;

    private PrintStream outputStream = null;


    @Override
    protected void processAlignments(JavaRDD<LocusWalkerContext> rdd, JavaSparkContext ctx) {
        rdd.map(intervalFunction(variants)).saveAsTextFile(outputFile);
    }

    private static Function<LocusWalkerContext, String> intervalFunction(List<FeatureInput<VariantContext>> variants) {
        return (Function<LocusWalkerContext, String>) context -> {
            AlignmentContext alignmentContext = context.getAlignmentContext();
            ReferenceContext referenceContext = context.getReferenceContext();
            FeatureContext featureContext = context.getFeatureContext();

            StringBuilder sb = new StringBuilder();

            // Get pileup and counts
            ReadPileup pileup = alignmentContext.getBasePileup();
            // print the locus and coverage
            sb.append(String.format("Current locus %s:%d (coverage=%s)\n", alignmentContext.getContig(),
                    alignmentContext.getPosition(), pileup.size()));
            // print the reference context if available
            if ( referenceContext.hasBackingDataSource() ) {
                sb.append("\tReference base(s): " + new String(referenceContext.getBases()));
                sb.append("\n");
            }
            // print the overlapping variants if there are some
            if(featureContext.hasBackingDataSource()) {
                List<VariantContext> vars = featureContext.getValues(variants);
                if(!vars.isEmpty()) {
                    sb.append("\tOverlapping variant(s):\n");
                    for (VariantContext variant : vars) {
                        sb.append(String.format("\t\t%s:%d-%d, Ref:%s, Alt(s):%s\n", variant.getContig(), variant.getStart(),
                                variant.getEnd(), variant.getReference(), variant.getAlternateAlleles()));
                    }
                }
            }

            return sb.toString();
        };
    }
}
