package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(summary="Align assembled contigs to the reference",
        oneLineSummary="Align assembled contigs to the reference",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class AlignAssembledContigsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    /**
     * Number of assemblies to process per task. We need to make sure that the work is not
     * over-partitioned since each worker needs to localize the BWA reference index and read it into
     * memory once per partition.
     */
    public static final int NUM_ASSEMBLIES_PER_PARTITION = 400;
    public static final int EXPECTED_CONTIGS_PER_ASSEMBLY = 15;

    @Argument(doc = "file for breakpoint alignment output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Argument(doc = "Input file of assembled contigs", shortName = "inputFile",
            fullName = "inputFile", optional = false)
    private String input;

    private static final Logger log = LogManager.getLogger(AlignAssembledContigsSpark.class);

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaPairRDD<String, ContigsCollection> breakpointIdsToContigsCollection = ContigsCollection.loadContigsCollectionKeyedByAssemblyId(ctx, input).cache();

        final long numInputAssemblies = breakpointIdsToContigsCollection.count();
        final int numPartitions = Math.max(ctx.defaultParallelism(), (int) Math.ceil((double) numInputAssemblies / (double) NUM_ASSEMBLIES_PER_PARTITION));

        final String referenceFileName = referenceArguments.getReferenceFileName();

        final JavaRDD<AlignmentRegion> allContigAlignments = breakpointIdsToContigsCollection.coalesce(numPartitions).mapPartitions(iter -> {
            try {
                try (final ContigAligner contigAligner = new ContigAligner(referenceFileName)) {
                    final List<AlignmentRegion> results = new ArrayList<>(NUM_ASSEMBLIES_PER_PARTITION * EXPECTED_CONTIGS_PER_ASSEMBLY);
                    iter.forEachRemaining(cc -> {
                        String breakpointId = cc._1;
                        final List<AlignmentRegion> contigAlignments = contigAligner.alignContigs(breakpointId, cc._2);
                        contigAlignments.forEach(results::add);
                    });
                    return results.iterator();
                }
            } catch (final IOException e) {
                throw new GATKException("Cannot run BWA-MEM", e);
            }

        });
        allContigAlignments.saveAsTextFile(output);
    }

}