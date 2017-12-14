package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.VariantWalkerContext;
import org.broadinstitute.hellbender.engine.spark.VariantWalkerSpark;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;

@CommandLineProgramProperties(summary = "Print variants from the input VCF", oneLineSummary = "Print variants from the input VCF", programGroup = SparkProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class PrintVariantsSpark extends VariantWalkerSpark {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String output;

    @Override
    protected void processVariants(JavaRDD<VariantWalkerContext> rdd, JavaSparkContext ctx) {
        try {
            VariantsSparkSink.writeVariants(ctx, output, rdd.map(VariantWalkerContext::getVariant), getHeaderForVariants());
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(output, "writing failed", e);
        }
    }
}
