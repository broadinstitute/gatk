package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.spark.VariantWalkerContext;
import org.broadinstitute.hellbender.engine.spark.VariantWalkerSpark;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;

/**
 * Print out variants from a VCF file.
 *
 * <h3>Input</h3>
 * <p>
 * A VCF file.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A new file containing the variants.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk PrintVariantsSpark \
 *    -V input.vcf.gz \
 *    -L chr1 \
 *    -O output.vcf.gz
 * </pre>
 *
 */

@CommandLineProgramProperties(
        summary = "Prints out variants from the input VCF file.",
        oneLineSummary = "Prints out variants from the input VCF.",
        programGroup = VariantManipulationProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class PrintVariantsSpark extends VariantWalkerSpark {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "Uri for the output file (a local file path).",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
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
