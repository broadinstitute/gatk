package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.spark.IntervalWalkerContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.engine.spark.IntervalWalkerSpark;

import java.io.PrintStream;
import java.util.List;

/**
 * Example/toy program that shows how to implement the IntervalWalker interface. Prints supplied intervals
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints intervals supplied via -L to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Print intervals with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleIntervalWalkerSpark extends IntervalWalkerSpark {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private OptionalVariantInputArgumentCollection optionalVariants = new OptionalVariantInputArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    private PrintStream outputStream = null;


    @Override
    protected void processIntervals(JavaRDD<IntervalWalkerContext> rdd, JavaSparkContext ctx) {
        rdd.map(intervalFunction(optionalVariants.variantFiles)).saveAsTextFile(outputFile);
    }

    private static Function<IntervalWalkerContext, String> intervalFunction(List<FeatureInput<VariantContext>> auxiliaryVariants) {
        return (Function<IntervalWalkerContext, String>) context -> {
            SimpleInterval interval = context.getInterval();
            ReadsContext readsContext = context.getReadsContext();
            ReferenceContext referenceContext = context.getReferenceContext();
            FeatureContext featureContext = context.getFeatureContext();

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
