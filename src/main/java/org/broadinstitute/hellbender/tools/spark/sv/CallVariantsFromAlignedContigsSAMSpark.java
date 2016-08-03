package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.IterableUtils;
import org.apache.commons.collections4.IteratorUtils;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AlignmentRegion;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.List;
import java.util.Optional;

import static org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.BreakpointAllele;

@CommandLineProgramProperties(summary="Parse SAM records and call SVs",
        oneLineSummary="Parse SAM records and call SVs",
        programGroup = SparkProgramGroup.class)
public class CallVariantsFromAlignedContigsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "URL of the output path", shortName = "outputPath",
            fullName = "outputPath", optional = false)
    private String outputPath;

    @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
            fullName = "fastaReference", optional = false)
    private String fastaReference;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return ReadFilterLibrary.MAPPED;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        Broadcast<ReferenceMultiSource> broadcastReference = ctx.broadcast(getReference());
        final JavaPairRDD<Tuple2<String, String>, Iterable<GATKRead>> alignmentsGroupedByName = getReads().mapToPair(r -> new Tuple2<>(new Tuple2<>(r.getName(), r.getName()), r)).groupByKey();
        final JavaPairRDD<Tuple2<String, String>, Tuple2<byte[], Iterable<AlignmentRegion>>> alignmentRegionsIterable = alignmentsGroupedByName.mapValues(CallVariantsFromAlignedContigsSAMSpark::convertToAlignmentRegions);
        final JavaPairRDD<Tuple2<String, String>, AssembledBreakpoint> alignedBreakpoints = alignmentRegionsIterable.flatMapValues(v -> CallVariantsFromAlignedContigsSpark.assembledBreakpointsFromAlignmentRegions(v._1, v._2));


        final JavaPairRDD<BreakpointAllele, Tuple2<Tuple2<String, String>, AssembledBreakpoint>> inversionBreakpoints =
                alignedBreakpoints
                        .mapToPair(CallVariantsFromAlignedContigsSpark::keyByBreakpointAllele)
                        .filter(CallVariantsFromAlignedContigsSpark::inversionBreakpointFilter);


        final JavaPairRDD<BreakpointAllele, Iterable<Tuple2<Tuple2<String,String>, AssembledBreakpoint>>> groupedBreakpoints = inversionBreakpoints.groupByKey();

        final JavaRDD<VariantContext> variantContexts = groupedBreakpoints.map(breakpoints -> CallVariantsFromAlignedContigsSpark.filterBreakpointsAndProduceVariants(breakpoints, broadcastReference)).cache();
        CallVariantsFromAlignedContigsSpark.writeVariants(fastaReference, logger, variantContexts, getAuthenticatedGCSOptions(), outputPath);
    }

    @VisibleForTesting
    protected static Tuple2<byte[], Iterable<AlignmentRegion>> convertToAlignmentRegions(final Iterable<GATKRead> reads) {
        final List<GATKRead> gatkReads = IterableUtils.toList(reads);
        final Optional<GATKRead> primaryRead = gatkReads.stream().filter(r -> !(r.isSecondaryAlignment() || r.isSupplementaryAlignment())).findFirst();
        if (! primaryRead.isPresent()) {
            throw new GATKException("no primary alignment for read " + gatkReads.get(0).getName());
        }

        final byte[] bases = primaryRead.get().getBases();
        if (primaryRead.get().isReverseStrand()) {
            SequenceUtil.reverseComplement(bases, 0, bases.length);
        }

        Iterable<AlignmentRegion> alignmentRegionIterable = IteratorUtils.asIterable(gatkReads.stream().filter(r -> !r.isSecondaryAlignment()).map(AlignmentRegion::new).iterator());
        return new Tuple2<>(bases, alignmentRegionIterable);
    }
}
