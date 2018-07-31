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
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionWalkerContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionWalkerSpark;

import java.io.Serializable;

/**
 * Example/toy program that shows how to implement the AssemblyRegionWalker interface.
 *
 * Prints out the bounds of each assembly region with and without padding, as well as the number of reads in each region.
 * Also prints overlapping reference bases / variants for each region, if provided.
 */
@CommandLineProgramProperties(
        summary = "Example AssemblyRegionWalker that prints out the bounds of each assembly region with and without padding, as well as the number of reads in each region",
        oneLineSummary = "Example AssemblyRegionWalker",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleAssemblyRegionWalkerSpark extends AssemblyRegionWalkerSpark {
    private static final long serialVersionUID = 1L;

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(fullName="knownVariants", shortName="knownVariants", doc="Known set of variants", optional=true)
    private FeatureInput<VariantContext> knownVariants;

    @Override
    protected int defaultReadShardSize() { return 5000; }

    @Override
    protected int defaultReadShardPadding() { return 100; }

    @Override
    protected int defaultMinAssemblyRegionSize() { return 50; }

    @Override
    protected int defaultMaxAssemblyRegionSize() { return 300; }

    @Override
    protected int defaultAssemblyRegionPadding() { return 100; }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() { return 50; }

    @Override
    protected double defaultActiveProbThreshold() { return 0.002; }

    @Override
    protected int defaultMaxProbPropagationDistance() { return 50; }

    @Override
    protected boolean includeReadsWithDeletionsInIsActivePileups() {
        return true;
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() {
        // This example AssemblyRegionEvaluator considers all loci to be active
        // (ie., it assigns an isActive probability of 1.0 to all loci):
        return new ARE();
    }

    static class ARE implements AssemblyRegionEvaluator, Serializable {
        private static final long serialVersionUID = 1L;
        @Override
        public ActivityProfileState isActive(AlignmentContext locusPileup, ReferenceContext referenceContext, FeatureContext featureContext) {
            return new ActivityProfileState(new SimpleInterval(locusPileup), 1.0);
        }
    }

    @Override
    protected void processAssemblyRegions(JavaRDD<AssemblyRegionWalkerContext> rdd, JavaSparkContext ctx) {
        rdd.map(assemblyFunction(knownVariants)).saveAsTextFile(outputFile);
    }

    private static Function<AssemblyRegionWalkerContext, String> assemblyFunction(FeatureInput<VariantContext> knownVariants) {
        return (Function<AssemblyRegionWalkerContext, String>) context -> {
            AssemblyRegion region = context.getAssemblyRegion();
            ReferenceContext referenceContext = context.getReferenceContext();
            FeatureContext featureContext = context.getFeatureContext();

            StringBuilder sb = new StringBuilder();
            sb.append(String.format("%s assembly region at %s (%s with padding), containing %d reads.\n\n",
                    region.isActive() ? "ACTIVE" : "INACTIVE", region.getSpan(), region.getExtendedSpan(), region.getReads().size()));
            sb.append(String.format("\tOverlapping reference bases: %s\n\n", new String(referenceContext.getBases())));

            if ( featureContext.hasBackingDataSource() ) {
                for ( final VariantContext variant : featureContext.getValues(knownVariants) ) {
                    sb.append(String.format("\tOverlapping variant at %s:%d-%d. Ref: %s Alt(s): %s\n\n",
                            variant.getContig(), variant.getStart(), variant.getEnd(), variant.getReference(), variant.getAlternateAlleles()));
                }
            }
            return sb.toString();
        };
    }

}
