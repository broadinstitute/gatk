package org.broadinstitute.hellbender.tools.spark.sv.playground;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark.annotateVariant;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark.inferType;


/**
 * This tool takes a SAM file containing alignments of single-ended long read
 * (be it long read sequencing, or contigs assembled from standard Illumina short reads),
 * searches for split alignments indicating the presence of structural variations,
 * and outputs an interpreted annotated single sample VCF.
 */
@CommandLineProgramProperties(summary="Parses a SAM file containing long reads alignments, and outputs an interpreted annotated single sample VCF.",
        oneLineSummary="Parses a long read SAM file, and outputs single sample VCF.",
        usageExample = "InternalDiscoverVariantsFromFilteredContigAlignmentsSAMSpark \\" +
                "-I /path/to/my/dir/longReads.sam -O /path/to/my/dir/structuralVariants.vcf \\" +
                "-R /path/to/my/reference/reference.2bit --fastaReference /path/to/my/reference/reference.fasta",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public class InternalDiscoverVariantsFromFilteredContigAlignmentsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(InternalDiscoverVariantsFromFilteredContigAlignmentsSAMSpark.class);

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
            discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "sam file for aligned contigs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutputFileName;

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

        // filter alignments and split the gaps
        final JavaRDD<AlignedContig> nameAndSeqAndFilteredAlignmentsAndGapSplit
                = InternalFilterLongReadAlignmentsSAMSpark.newWayOfFiltering(getReads(), getHeaderForReads(), localLogger);

        // convert to ChimericAlignment, similar to ChimericAlignment.parseOneContig(final AlignedContig, final int), except the head/tail filtering
        final JavaPairRDD<byte[], List<ChimericAlignment>> chimericAlignments =
                nameAndSeqAndFilteredAlignmentsAndGapSplit
                        .filter(contig -> contig.alignmentIntervals.size() > 1) // remove contigs who after filtering has only one alignment
                        .mapToPair(InternalDiscoverVariantsFromFilteredContigAlignmentsSAMSpark::convertAlignmentIntervalToChimericAlignment);

        inferSvAndWriteVCF(chimericAlignments, vcfOutputFileName, ctx.broadcast(getReference()),
                discoverStageArgs.fastaReference, getAuthenticatedGCSOptions(), localLogger);
    }

    /**
     * Very similar to {@link ChimericAlignment#parseOneContig(AlignedContig, int)}, except the head/tail filtering.
     */
    private static Tuple2<byte[], List<ChimericAlignment>> convertAlignmentIntervalToChimericAlignment
            (final AlignedContig contig) {

        final List<AlignmentInterval> alignmentIntervals = contig.alignmentIntervals;
        final Iterator<AlignmentInterval> iterator = alignmentIntervals.iterator();
        AlignmentInterval current = iterator.next();
        final List<ChimericAlignment> results = new ArrayList<>(alignmentIntervals.size() - 1);
        final List<String> insertionMappings = new ArrayList<>();
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            if (nextAlignmentMayBeNovelInsertion(current, next, DEFAULT_MIN_ALIGNMENT_LENGTH)) {
                if (iterator.hasNext()) {
                    insertionMappings.add(next.toPackedString());
                    continue;
                } else {
                    break;
                }
            }
            final boolean isNotSimpleTranslocation = isNotSimpleTranslocation(current, next,
                    determineStrandSwitch(current, next), involvesRefPositionSwitch(current, next));
            if (isNotSimpleTranslocation)
                results.add(new ChimericAlignment(current, next, insertionMappings, contig.contigName));
        }
        return new Tuple2<>(contig.contigSequence,results);
    }

    // usual business as in DiscoverVariantsFromContigAlignmentsSAMSpark#discoverVariantsAndWriteVCF()
    private static void inferSvAndWriteVCF(final JavaPairRDD<byte[], List<ChimericAlignment>> chimericAlignments, final String vcfOutputFileName,
                                           final Broadcast<ReferenceMultiSource> broadcastReference, final String fastaReference,
                                           final GCSOptions options, final Logger toolLogger) {

        final JavaRDD<VariantContext> annotatedVariants =
                chimericAlignments
                        .flatMapToPair(DiscoverVariantsFromContigAlignmentsSAMSpark::discoverNovelAdjacencyFromChimericAlignments)
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> inferType(noveltyAndEvidence._1, noveltyAndEvidence._2))
                        .map(noveltyTypeAndEvidence -> annotateVariant(noveltyTypeAndEvidence._1,
                                noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference));

        SVVCFWriter.writeVCF(options, vcfOutputFileName, fastaReference, annotatedVariants, toolLogger);
    }
}
