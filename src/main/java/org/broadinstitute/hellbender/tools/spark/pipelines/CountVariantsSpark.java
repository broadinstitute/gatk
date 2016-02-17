package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

@CommandLineProgramProperties(summary = "Counts variants in the input VCF",
        oneLineSummary = "CountVariants on Spark",
        programGroup = SparkProgramGroup.class)
public final class CountVariantsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the input file: a local file path",
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            optional = true)
    public String input;


    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final VariantsSparkSource vss = new VariantsSparkSource(ctx);
        final JavaRDD<VariantContext> variants = vss.getParallelVariantContexts(input);

        System.out.println( variants.first() );

        final long count = variants.count();
        System.out.println(count);

        if (out != null){
            final File file = new File(out);
            try(final OutputStream outputStream = BucketUtils.createNonGCSFile(file.getPath());
                final PrintStream ps = new PrintStream(outputStream)) {
                ps.print(count);
            } catch(final IOException e){
                throw new UserException.CouldNotCreateOutputFile(file, e);
            }
        }

    }
}
