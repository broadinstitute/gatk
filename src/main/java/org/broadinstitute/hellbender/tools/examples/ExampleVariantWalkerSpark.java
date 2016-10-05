package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.engine.spark.VariantWalkerSpark;
import scala.Tuple4;

import java.io.PrintStream;

/**
 * Example/toy program that shows how to implement the VariantWalker interface. Prints supplied variants
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints variants supplied to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Example tool that prints variants with optional contextual data",
        programGroup = VariantProgramGroup.class
)
public final class ExampleVariantWalkerSpark extends VariantWalkerSpark {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    private PrintStream outputStream = null;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        getVariants(ctx).map(variantFunction(auxiliaryVariants)).saveAsTextFile(outputFile);
    }

    private static Function<Tuple4<VariantContext, ReadsContext, ReferenceContext, FeatureContext>, String> variantFunction(FeatureInput<VariantContext> auxiliaryVariants) {
        return (Function<Tuple4<VariantContext, ReadsContext, ReferenceContext, FeatureContext>, String>) t -> {
            VariantContext variant = t._1();
            ReadsContext readsContext = t._2();
            ReferenceContext referenceContext = t._3();
            FeatureContext featureContext = t._4();

            StringBuilder sb = new StringBuilder();
            sb.append(String.format("Current variant: " + variant));
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
