package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.spark.VariantWalkerContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.engine.spark.VariantWalkerSpark;

import java.io.PrintStream;

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
public final class ExampleVariantWalkerSpark extends VariantWalkerSpark {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    @Override
    protected void processVariants(JavaRDD<VariantWalkerContext> rdd, JavaSparkContext ctx) {
        rdd.map(variantFunction(auxiliaryVariants)).saveAsTextFile(outputFile);
    }

    private static Function<VariantWalkerContext, String> variantFunction(FeatureInput<VariantContext> auxiliaryVariants) {
        return (Function<VariantWalkerContext, String>) context -> {
            VariantContext variant = context.getVariant();
            ReadsContext readsContext = context.getReadsContext();
            ReferenceContext referenceContext = context.getReferenceContext();
            FeatureContext featureContext = context.getFeatureContext();

            StringBuilder sb = new StringBuilder();
            sb.append("Current variant: " + variant);
            sb.append("\n");

            if ( referenceContext.hasBackingDataSource() ) {
                sb.append(String.format("\tOverlapping reference bases: %s\n\n", new String(referenceContext.getBases())));
            }

            if ( readsContext.hasBackingDataSource() ) {
                for ( final GATKRead read : readsContext) {
                    sb.append(String.format("\tOverlapping read at %s:%d-%d\n", read.getContig(), read.getStart(), read.getEnd()));
                }
                sb.append("\n");
            }

            if ( featureContext.hasBackingDataSource() ) {
                for ( final VariantContext variant1 : featureContext.getValues(auxiliaryVariants) ) {
                    sb.append(String.format("\tOverlapping variant at %s:%d-%d. Ref: %s Alt(s): %s\n",
                            variant1.getContig(), variant1.getStart(), variant1.getEnd(), variant1.getReference(), variant1.getAlternateAlleles()));
                }
                sb.append("\n");
            }

            return sb.toString();
        };
    }
}
