package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.IntervalProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.engine.spark.IntervalWalkerSpark;
import scala.Tuple4;

import java.io.PrintStream;
import java.util.List;

/**
 * Example/toy program that shows how to implement the IntervalWalker interface. Prints supplied intervals
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints intervals supplied via -L to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Print intervals with optional contextual data",
        programGroup = IntervalProgramGroup.class
)
public final class ExampleIntervalWalkerSpark extends IntervalWalkerSpark {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private OptionalVariantInputArgumentCollection optionalVariants = new OptionalVariantInputArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    private PrintStream outputStream = null;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        getIntervals(ctx).map(intervalFunction(optionalVariants.variantFiles)).saveAsTextFile(outputFile);
    }

    private static Function<Tuple4<SimpleInterval, ReadsContext, ReferenceContext, FeatureContext>, String> intervalFunction(List<FeatureInput<VariantContext>> auxiliaryVariants) {
        return (Function<Tuple4<SimpleInterval, ReadsContext, ReferenceContext, FeatureContext>, String>) t -> {
            SimpleInterval interval = t._1();
            ReadsContext readsContext = t._2();
            ReferenceContext referenceContext = t._3();
            FeatureContext featureContext = t._4();

            StringBuilder sb = new StringBuilder();
            sb.append(String.format("Current interval: " + interval));
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
                for ( final VariantContext variant : featureContext.getValues(auxiliaryVariants) ) {
                    sb.append(String.format("\tOverlapping variant at %s:%d-%d. Ref: %s Alt(s): %s\n",
                            variant.getContig(), variant.getStart(), variant.getEnd(), variant.getReference(), variant.getAlternateAlleles()));
                }
                sb.append("\n");
            }

            return sb.toString();
        };
    }
}
