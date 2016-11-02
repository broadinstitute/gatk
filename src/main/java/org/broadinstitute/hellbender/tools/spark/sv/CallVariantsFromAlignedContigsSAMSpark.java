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
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.List;
import java.util.Optional;

import static org.broadinstitute.hellbender.tools.spark.sv.CallVariantsFromAlignedContigsSpark.*;

/**
 * This tool takes a SAM file containing the alignments of assembled contigs or long reads to the reference
 * and searches it for split alignments indicating the presence of structural variations. To do so the tool parses
 * primary and supplementary alignments; secondary alignments are ignored. To be considered valid evidence of an SV,
 * two alignments from the same contig must have mapping quality 60, and both alignments must have length greater than
 * or equal to minAlignmentLength.
 */
@CommandLineProgramProperties(summary="Parse a SAM file containing contigs or long reads aligned to the reference, and call SVs",
        oneLineSummary="Parse a SAM file containing contigs or long reads aligned to the reference, and call SVs",
        programGroup = SparkProgramGroup.class)
public final class CallVariantsFromAlignedContigsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "URL of the output path", shortName = "outputPath",
            fullName = "outputPath", optional = false)
    private String outputPath;

    @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
            fullName = "fastaReference", optional = false)
    private String fastaReference;

    @Argument(doc = "Minimum flanking alignment length", shortName = "minAlignLength",
            fullName = "minAlignLength", optional = true)
    private Integer minAlignLength = SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH;

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
        final Broadcast<ReferenceMultiSource> broadcastReference = ctx.broadcast(getReference());
        final JavaPairRDD<Tuple2<String, String>, Iterable<GATKRead>> alignmentsGroupedByName = getReads().mapToPair(r -> new Tuple2<>(new Tuple2<>(r.getName(), r.getName()), r)).groupByKey();
        final JavaPairRDD<Tuple2<String, String>, Tuple2<Iterable<AlignmentRegion>, byte[]>> alignmentRegionsIterable = alignmentsGroupedByName.mapValues(CallVariantsFromAlignedContigsSAMSpark::convertToAlignmentRegions);

        final Integer minAlignLengthFinal = minAlignLength;

        final JavaRDD<VariantContext> variants = callVariantsFromAlignmentRegions(broadcastReference, alignmentRegionsIterable);
        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputPath, SVConstants.INVERSIONS_OUTPUT_VCF, fastaReference, variants);
    }

    @VisibleForTesting
    static Tuple2<Iterable<AlignmentRegion>, byte[]> convertToAlignmentRegions(final Iterable<GATKRead> reads) {
        final List<GATKRead> gatkReads = IterableUtils.toList(reads);
        final Optional<GATKRead> primaryRead = gatkReads.stream().filter(r -> !(r.isSecondaryAlignment() || r.isSupplementaryAlignment())).findFirst();
        if (! primaryRead.isPresent()) {
            throw new GATKException("no primary alignment for read " + gatkReads.get(0).getName());
        }

        final byte[] bases = primaryRead.get().getBases();
        if (primaryRead.get().isReverseStrand()) {
            SequenceUtil.reverseComplement(bases, 0, bases.length);
        }

        Iterable<AlignmentRegion> alignmentRegionIterable = IteratorUtils.asIterable(gatkReads.stream()
                .filter(r -> !r.isSecondaryAlignment()).map(AlignmentRegion::new).iterator());
        return new Tuple2<>(alignmentRegionIterable, bases);
    }
}
