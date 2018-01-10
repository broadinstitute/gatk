package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(summary = "Counts variants in the input VCF",
        oneLineSummary = "CountVariants on Spark",
        programGroup = VariantEvaluationProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class CountVariantsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the input file: a local file path",
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            optional = false)
    public String input;


    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    /**
     * Return the list of GATKCommandLinePluginDescriptor objects to be used for this CLP.
     * GATKSparkTool returns the GATKReadFilterPluginDescriptor, but we don't want that
     * for Variant tools.
     */
    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        return new ArrayList<>();
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final VariantsSparkSource vss = new VariantsSparkSource(ctx);
        final JavaRDD<VariantContext> variants = vss.getParallelVariantContexts(input, getIntervals());

        final long count = variants.count();
        System.out.println(count);

        if( out != null) {
            try (final PrintStream ps = new PrintStream(BucketUtils.createFile(out))) {
                ps.print(count);
            }
        }
    }
}
