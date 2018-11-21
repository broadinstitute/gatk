package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.spark.ReadWalkerContext;
import org.broadinstitute.hellbender.engine.spark.ReadWalkerSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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
public final class ExampleReadWalkerWithVariantsSpark extends ReadWalkerSpark {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files", optional = true)
    private List<FeatureInput<VariantContext>> variants;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    @Argument(fullName = "groupVariantsBySource", doc = "If true, group overlapping variants by their source when outputting them", optional = true)
    private boolean groupVariantsBySource = false;

    @Override
    protected void processReads(JavaRDD<ReadWalkerContext> rdd, JavaSparkContext ctx) {
        rdd.map(readFunction(variants, groupVariantsBySource)).saveAsTextFile(outputFile);
    }

    private static Function<ReadWalkerContext, String> readFunction(List<FeatureInput<VariantContext>> variants, boolean groupVariantsBySource) {
        return (Function<ReadWalkerContext, String>) context -> {
            GATKRead read = context.getRead();
            ReferenceContext referenceContext = context.getReferenceContext();
            FeatureContext featureContext = context.getFeatureContext();

            StringBuilder sb = new StringBuilder();

            sb.append(String.format("Read at %s:%d-%d:\n%s\n", read.getContig(), read.getStart(), read.getEnd(), read.getBasesString()));

            if ( groupVariantsBySource ) {
                // We can keep the variants from each source separate by passing in the FeatureInputs
                // individually to featureContext.getValues()
                for ( FeatureInput<VariantContext> featureSource : variants ) {
                    sb.append("From source " + featureSource.getName());
                    sb.append("\n");

                    for ( VariantContext variant : featureContext.getValues(featureSource) ) {
                        sb.append(String.format("\t"));
                        printOverlappingVariant(sb, variant);
                    }
                }
            }
            else {
                // Passing in all FeatureInputs at once to featureContext.getValues() lets us get
                // all overlapping variants without regard to source
                for ( VariantContext variant : featureContext.getValues(variants) ) {
                    printOverlappingVariant(sb, variant);
                }
            }

            sb.append("\n");

            return sb.toString();
        };
    }

    private static void printOverlappingVariant(StringBuilder sb, final VariantContext variant) {
        sb.append(String.format("Overlapping variant at %s:%d-%d: Ref: %s Alt(s): %s\n",
                variant.getContig(), variant.getStart(), variant.getEnd(),
                variant.getReference(), variant.getAlternateAlleles()));
    }
}
