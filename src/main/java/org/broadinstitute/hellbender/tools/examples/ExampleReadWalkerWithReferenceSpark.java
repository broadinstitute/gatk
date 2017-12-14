package org.broadinstitute.hellbender.tools.examples;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.spark.ReadWalkerContext;
import org.broadinstitute.hellbender.engine.spark.ReadWalkerSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Example/toy program that prints reads from the provided file or files with corresponding reference bases
 * (if a reference is provided). Intended to show how to implement the ReadWalker interface.
 */
@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) with corresponding reference bases (if a reference is provided) to the specified output file (or STDOUT if none specified)",
        oneLineSummary = "Print reads with reference context",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class ExampleReadWalkerWithReferenceSpark extends ReadWalkerSpark {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    protected void processReads(JavaRDD<ReadWalkerContext> rdd, JavaSparkContext ctx) {
        rdd.map(readFunction()).saveAsTextFile(outputFile);
    }

    private Function<ReadWalkerContext, String> readFunction() {
        return (Function<ReadWalkerContext, String>) context -> {
            GATKRead read = context.getRead();
            ReferenceContext referenceContext = context.getReferenceContext();

            StringBuilder sb = new StringBuilder();
            sb.append(String.format("Read at %s:%d-%d:\n%s\n", read.getContig(), read.getStart(), read.getEnd(), read.getBasesString()));
            if ( referenceContext.hasBackingDataSource() )
                sb.append("Reference Context:\n" + new String(referenceContext.getBases()) + "\n");
            sb.append("\n");

            return sb.toString();
        };
    }
}
