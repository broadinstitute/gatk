package org.broadinstitute.hellbender.tools.spark.sv.sga;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.lang3.StringUtils;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedAssembly;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexSingleton;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary="Align assembled contigs to the reference",
        oneLineSummary="Align assembled contigs to the reference",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class AlignAssembledContigsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    public static final String MAPPED_CONTIG_ALIGNMENT_INTERVAL_STRING_REP_FIELD_SEPARATOR = "%";
    public static final String UNMAPPED_CONTIG_STRING_REP = "unmapped";

    /**
     * Number of assemblies to process per task. We need to make sure that the work is not
     * over-partitioned since each worker needs to localize the BWA reference index and read it into
     * memory once per partition.
     */
    private static final int NUM_ASSEMBLIES_PER_PARTITION = 400;
    private static final int EXPECTED_CONTIGS_PER_ASSEMBLY = 15;

    @Argument(doc = "Input directory of assembled contigs", shortName = "inputAssemblyDir", fullName = "inputAssemblyDir")
    private String input;

    @Argument(doc = "file for breakpoint alignment output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String output;

    @Argument(doc = "the bwa mem index image file name that you've distributed to each executor", fullName = "bwamemIndexImage")
    private String indexImageFile;

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaPairRDD<String, ContigsCollection> breakpointIdsToContigsCollection = ContigsCollection.loadContigsCollectionKeyedByAssemblyId(ctx, input).cache();

        final long numInputAssemblies = breakpointIdsToContigsCollection.count();
        final int numPartitions = Math.max(ctx.defaultParallelism(), (int) Math.ceil((double) numInputAssemblies / (double) NUM_ASSEMBLIES_PER_PARTITION));
        final String indexImageFile = this.indexImageFile;

        final JavaRDD<AlignedAssembly> allContigAlignments =
                breakpointIdsToContigsCollection
                        .coalesce(numPartitions)
                        .mapPartitions(iter -> {
                            final ContigAligner contigAligner = new ContigAligner(indexImageFile);
                            final List<AlignedAssembly> alignedAssembliesInThisPartition = new ArrayList<>(NUM_ASSEMBLIES_PER_PARTITION);
                            iter.forEachRemaining(pair -> alignedAssembliesInThisPartition.add(contigAligner.alignContigs(Integer.valueOf(pair._1), pair._2)));
                            return alignedAssembliesInThisPartition.iterator();
                        });

        allContigAlignments.flatMap(AlignAssembledContigsSpark::formatAlignedAssemblyAsText).saveAsTextFile(output);

        BwaMemIndexSingleton.closeAllDistributedInstances(ctx);
    }

    /**
     * Format input aligned assembly as list of strings where each entry is for one aligned contig.
     * And each aligned contig is formatted as:
     * <pre>
     * contigName ({@link AlignedAssemblyOrExcuse#formatContigName(int, int)}) + TAB + {ALIGNMENT_INTERVAL}
     * </pre>
     * where {ALIGNMENT_INTERVAL} is a list of formatted {@link AlignedAssembly.AlignmentInterval}'s separated by TAB's,
     * and each {@link AlignedAssembly.AlignmentInterval} is formatted as:
     * <pre>
     * startInAssembledContig + "-" + endInAssembledContig + {@link #MAPPED_CONTIG_ALIGNMENT_INTERVAL_STRING_REP_FIELD_SEPARATOR} +
     * referenceInterval({@link org.broadinstitute.hellbender.utils.SimpleInterval#toString()} + {@link #MAPPED_CONTIG_ALIGNMENT_INTERVAL_STRING_REP_FIELD_SEPARATOR} +
     * cigarAlong5to3DirectionOfContig + {@link #MAPPED_CONTIG_ALIGNMENT_INTERVAL_STRING_REP_FIELD_SEPARATOR} +
     * "+/-" (depending on strandedness) + {@link #MAPPED_CONTIG_ALIGNMENT_INTERVAL_STRING_REP_FIELD_SEPARATOR} +
     * mapQual + {@link #MAPPED_CONTIG_ALIGNMENT_INTERVAL_STRING_REP_FIELD_SEPARATOR} + mismatch
     * </pre>
     */
    @VisibleForTesting
    public static Iterator<String> formatAlignedAssemblyAsText(final AlignedAssembly alignedAssembly) {

        return alignedAssembly.alignedContigs.stream()
                .map(alignedContig -> {
                    final String mappingInfo;
                    if (alignedContig.alignmentIntervals.isEmpty()) {
                        mappingInfo = UNMAPPED_CONTIG_STRING_REP;
                    } else {
                        final List<String> intervals =
                                alignedContig.alignmentIntervals.stream()
                                        .map(alignmentInterval ->
                                                StringUtils.join(Arrays.asList(String.valueOf(alignmentInterval.startInAssembledContig) + "-" + String.valueOf(alignmentInterval.endInAssembledContig), encodeSimpleIntervalAsString(alignmentInterval.referenceInterval),TextCigarCodec.encode(alignmentInterval.cigarAlong5to3DirectionOfContig), (alignmentInterval.forwardStrand ? "+" : "-"), alignmentInterval.mapQual, alignmentInterval.mismatches), MAPPED_CONTIG_ALIGNMENT_INTERVAL_STRING_REP_FIELD_SEPARATOR))
                                        .collect(Collectors.toList());
                        mappingInfo = StringUtils.join(intervals,"\t");
                    }
                    return alignedContig.contigName + "\t" + mappingInfo;
                })
                .collect(Collectors.toList()).iterator();
    }

    /**
     * Preferably, we should use {@link SimpleInterval#toString()}, but HG38 contig names such as
     * <pre>
     *     ##contig=&lt;ID=HLA-DRB1*16:02:01,length=11005,assembly=38&gt;
     * </pre>
     * forbids such move.
     */
    public static String encodeSimpleIntervalAsString(final SimpleInterval simpleInterval) {
        return "CTG=" + simpleInterval.getContig() + "START=" + String.valueOf(simpleInterval.getStart()) + "END=" + String.valueOf(simpleInterval.getEnd());
    }

    @VisibleForTesting
    public static SimpleInterval decodeStringAsSimpleInterval(final String text) {
        try {
            if (!text.startsWith("CTG=")) throw new IllegalArgumentException(text);
            int stop1 = text.indexOf("START=");
            final String ctg = text.substring(4, stop1);
            int stop2 = text.indexOf("END=", stop1);
            final int start = Integer.valueOf(text.substring(stop1+6, stop2));
            final int end = Integer.valueOf(text.substring(stop2+4, text.length()));
            return new SimpleInterval(ctg, start, end);
        } catch (final Exception ex) {
            throw new GATKException("unexpected format in supposedly formatted text for a SimpleInterval: " + text, ex);
        }
    }

}
