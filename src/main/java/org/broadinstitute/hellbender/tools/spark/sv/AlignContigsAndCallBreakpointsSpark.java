package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AssembledBreakpoint;
import org.broadinstitute.hellbender.tools.spark.sv.RunSGAViaProcessBuilderOnSpark.ContigsCollection;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@CommandLineProgramProperties(summary="Align assembled contigs to the reference and call breakpoints from them.",
        oneLineSummary="Align contigs and call breakpoints",
        programGroup = SparkProgramGroup.class)
public class AlignContigsAndCallBreakpointsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    /**
     * Number of assemblies to process per task. We need to make sure that the work is not
     * over-partitioned since each worker needs to localize the BWA reference index and read it into
     * memory once per partition.
     */
    public static final int NUM_ASSEMBLIES_PER_PARTITION = 400;
    public static final int EXPECTED_BREAKPOINTS_PER_ASSEMBLY = 10;

    @Argument(doc = "file for assembled breakpoint output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Argument(doc = "Input file of assembled contigs", shortName = "inputFile",
            fullName = "inputFile", optional = false)
    private String input;

    private static final Logger log = LogManager.getLogger(AlignContigsAndCallBreakpointsSpark.class);

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<String> inputAssemblies = ctx.textFile(input).cache();

        final long numInputPartitions = inputAssemblies.count();

        final int numPartitions = Math.max(ctx.defaultParallelism (), (int) Math.ceil((double) numInputPartitions / (double) NUM_ASSEMBLIES_PER_PARTITION));
        final JavaPairRDD<String, String> contigCollectionByBreakpointId =
                inputAssemblies
                        .flatMapToPair(RunSGAViaProcessBuilderOnSpark::splitAssemblyLine)
                        .coalesce(numPartitions);

        final JavaPairRDD<String, ContigsCollection> breakpointIdsToContigsCollection =
                contigCollectionByBreakpointId.mapValues(ContigsCollection::fromPackedFasta);
        final String referenceFileName = referenceArguments.getReferenceFileName();

        final JavaPairRDD<String, AssembledBreakpoint> assembledBreakpoints = breakpointIdsToContigsCollection.mapPartitionsToPair(iter -> {
            try {
                try (final ContigAligner contigAligner = new ContigAligner(referenceFileName)) {
                    final List<Tuple2<String, AssembledBreakpoint>> results = new ArrayList<>(NUM_ASSEMBLIES_PER_PARTITION * EXPECTED_BREAKPOINTS_PER_ASSEMBLY);
                    iter.forEachRemaining(cc -> {
                        final List<AssembledBreakpoint> breakpointAssemblies = contigAligner.alignContigs(cc._2);
                        breakpointAssemblies.forEach(b -> results.add(new Tuple2<>(cc._1, b)));
                    });
                    return results;
                }
            } catch (final IOException e) {
                throw new GATKException("Cannot run BWA-MEM", e);
            }

        });
        assembledBreakpoints.saveAsTextFile(output);
    }

    /**
     * input format is comma-separated BreakpointId, String representation of an AssembledBreakpoint, enclosed in parens
     * @param alignedAssemblyContigLine An input line with a breakpoint ID and string representation of an AssembledBreakpoint
     * @return A tuple with the breakpoint ID and string representation of an AssembledBreakpoint, or an empty iterator if the line did not have two comma-separated values
     */
    static Iterable<Tuple2<String, AssembledBreakpoint>> parseAlignedAssembledContigLine(final String alignedAssemblyContigLine) {
        final String[] split = alignedAssemblyContigLine.replace("(","").replace(")","").split(",");
        if (split.length < 2) {
            log.info("No aligned breakpoints for line " + alignedAssemblyContigLine);
            return Collections.emptySet();
        }
        return Collections.singleton(new Tuple2<>(split[0], AssembledBreakpoint.fromFields(split[1].split("\t", -1))));
    }

}