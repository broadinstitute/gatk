package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AlignmentRegion;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AssembledBreakpoint;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.BreakpointAllele;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary="Filter breakpoint alignments and call variants.",
        oneLineSummary="Filter breakpoint alignments and call variants",
        programGroup = SparkProgramGroup.class)
public class CallVariantsFromAssembledBreakpointsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "Called variants output file", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Argument(doc = "Input file of aligned assembled contigs", shortName = "inputFile",
            fullName = "inputFile", optional = false)
    private String input;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaPairRDD<String, AssembledBreakpoint> inputAlignedBreakpoints = ctx.textFile(input).flatMapToPair(AlignContigsAndCallBreakpointsSpark::parseAlignedAssembledContigLine).cache();

        final JavaPairRDD<BreakpointAllele, Tuple2<String, AssembledBreakpoint>> assembled3To5BreakpointsKeyedByPosition =
                inputAlignedBreakpoints
                        .filter(CallVariantsFromAssembledBreakpointsSpark::inversion3To5BreakpointFilter)
                        .mapToPair(CallVariantsFromAssembledBreakpointsSpark::keyByBreakpointAllele);

        final JavaPairRDD<BreakpointAllele, Iterable<Tuple2<String, AssembledBreakpoint>>> groupedBreakpoints = assembled3To5BreakpointsKeyedByPosition.groupByKey();

        final JavaRDD<VariantContext> variantContexts = groupedBreakpoints.map(CallVariantsFromAssembledBreakpointsSpark::filterBreakpointsAndProduceVariants);

        variantContexts.saveAsTextFile(output);

    }

    @VisibleForTesting
    private static VariantContext filterBreakpointsAndProduceVariants(final Tuple2<BreakpointAllele, Iterable<Tuple2<String, AssembledBreakpoint>>> assembledBreakpointsPerAllele) throws IOException {
        int numAssembledBreakpoints = 0;
        int highMqMappings = 0;
        int midMqMappings = 0;
        int lowMqMappings = 0;
        final List<Integer> mqs = new ArrayList<>(10);
        final List<Integer> alignLengths = new ArrayList<>(10);
        final BreakpointAllele breakpointAllele = assembledBreakpointsPerAllele._1;
        final Iterable<Tuple2<String, AssembledBreakpoint>> assembledBreakpoints = assembledBreakpointsPerAllele._2;
        for (final Tuple2<String, AssembledBreakpoint> assembledBreakpointPair : assembledBreakpoints) {
            final AssembledBreakpoint assembledBreakpoint = assembledBreakpointPair._2;
            numAssembledBreakpoints = numAssembledBreakpoints + 1;
            final int assembledBreakpointMapq = Math.min(assembledBreakpoint.region1.mqual, assembledBreakpoint.region2.mqual);
            if (assembledBreakpointMapq == 60) {
                highMqMappings = highMqMappings + 1;
            } else if (assembledBreakpointMapq > 0) {
                midMqMappings = midMqMappings + 1;
            } else {
                lowMqMappings = lowMqMappings + 1;
            }
            mqs.add(assembledBreakpointMapq);
            final int assembledBreakpointAlignmentLength = Math.min(assembledBreakpoint.region1.referenceInterval.size(), assembledBreakpoint.region2.referenceInterval.size());
            alignLengths.add(assembledBreakpointAlignmentLength);

        }

        return createVariant(numAssembledBreakpoints, highMqMappings, mqs, alignLengths, breakpointAllele);
    }

    @VisibleForTesting
    private static VariantContext createVariant(final int numAssembledBreakpoints, final int highMqMappings, final List<Integer> mqs, final List<Integer> alignLengths, final BreakpointAllele breakpointAllele) throws IOException {
        final Allele refAllele = Allele.create(refBase, true);
        final Allele altAllele = Allele.create("<INV>");
        final List<Allele> vcAlleles = new ArrayList<>(2);
        vcAlleles.add(refAllele);
        vcAlleles.add(altAllele);

        final VariantContextBuilder vcBuilder = new VariantContextBuilder()
            .chr(breakpointAllele.leftAlignedLeftBreakpoint.getContig())
            .start(breakpointAllele.leftAlignedLeftBreakpoint.getStart())
            .stop(breakpointAllele.leftAlignedRightBreakpoint.getStart())
            .alleles(vcAlleles)
            .attribute("SVTYPE", "INV")
            .attribute("GATKSVTYPE", "INV3")
            .attribute("TOTAL_MAPPINGS", numAssembledBreakpoints)
            .attribute("HQMAPPINGS", highMqMappings)
            .attribute("MQS", mqs.stream().map(String::valueOf).collect(Collectors.joining(",")))
            .attribute("ALIGN_LENGTHS", alignLengths.stream().map(String::valueOf).collect(Collectors.joining(",")));
        return vcBuilder.make();
    }

    private static boolean inversion3To5BreakpointFilter(final Tuple2<String, AssembledBreakpoint> assembledBreakpoint) {
        final AlignmentRegion region1 = assembledBreakpoint._2.region1;
        final AlignmentRegion region2 = assembledBreakpoint._2.region2;
        return region1.forwardStrand && ! region2.forwardStrand;
    }

    private static Tuple2<BreakpointAllele, Tuple2<String, AssembledBreakpoint>> keyByBreakpointAllele(final Tuple2<String, AssembledBreakpoint> breakpointIdAndAssembledBreakpoint) {
        final BreakpointAllele breakpointAllele = breakpointIdAndAssembledBreakpoint._2.getBreakpointAllele();
        return new Tuple2<>(breakpointAllele, new Tuple2<>(breakpointIdAndAssembledBreakpoint._1, breakpointIdAndAssembledBreakpoint._2));
    }
}
