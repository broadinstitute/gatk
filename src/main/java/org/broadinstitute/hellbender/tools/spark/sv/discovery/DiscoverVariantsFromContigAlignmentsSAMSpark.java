package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;


/**
 * This tool takes a SAM file containing the alignments of assembled contigs or long reads to the reference
 * and searches it for split alignments indicating the presence of structural variations. To do so the tool parses
 * primary and supplementary alignments; secondary alignments are ignored. To be considered valid evidence of an SV,
 * two alignments from the same contig must have mapping quality 60, and both alignments must have length greater than
 * or equal to minAlignmentLength.
 */
@CommandLineProgramProperties(summary="Parse a SAM file containing contigs or long reads aligned to the reference, and call SVs",
        oneLineSummary="Parse a SAM file containing contigs or long reads aligned to the reference, and call SVs",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class DiscoverVariantsFromContigAlignmentsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(DiscoverVariantsFromContigAlignmentsSAMSpark.class);


    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
            discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "output path for discovery (non-genotyped) VCF", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
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

        final SAMFileHeader headerForReads = getHeaderForReads();

        final JavaRDD<AlignedContig> parsedContigAlignments
                = new SAMFormattedContigAlignmentParser(getReads(), headerForReads, true, localLogger)
                .getAlignedContigs();

        discoverVariantsAndWriteVCF(parsedContigAlignments,
                discoverStageArgs,
                ctx.broadcast(getReference()),
                vcfOutputFileName,
                localLogger,
                headerForReads.getSequenceDictionary());
    }

    public static final class SAMFormattedContigAlignmentParser extends AlignedContigGenerator implements Serializable {
        private static final long serialVersionUID = 1L;

        private final JavaRDD<GATKRead> unfilteredContigAlignments;
        private final SAMFileHeader header;
        private final boolean splitGapped;
        private final Logger toolLogger;

        public SAMFormattedContigAlignmentParser(final JavaRDD<GATKRead> unfilteredContigAlignments,
                                                 final SAMFileHeader header, final boolean splitGapped,
                                                 final Logger toolLogger) {
            this.unfilteredContigAlignments = unfilteredContigAlignments;
            this.header = header;
            this.toolLogger = toolLogger;
            this.splitGapped = splitGapped;
        }

        @Override
        public JavaRDD<AlignedContig> getAlignedContigs() {
            return unfilteredContigAlignments
                    .filter(r -> !r.isSecondaryAlignment())
                    .groupBy(GATKRead::getName)
                    .map(Tuple2::_2)
                    .map(gatkReads ->
                            parseReadsAndOptionallySplitGappedAlignments(
                                    Utils.stream(gatkReads).map(r->r.convertToSAMRecord(header)).collect(Collectors.toList()),
                                    GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, splitGapped, toolLogger));
        }

        /**
         * Iterates through the input {@code noSecondaryAlignments}, which are assumed to contain no secondary alignment (i.e. records with "XA" tag),
         * converts to custom {@link AlignmentInterval} format and
         * split the records when the gap in the alignment reaches the specified {@code sensitivity}.
         * The size of the returned iterable of {@link AlignmentInterval}'s is guaranteed to be no lower than that of the input iterable.
         */
        @VisibleForTesting
        public static AlignedContig parseReadsAndOptionallySplitGappedAlignments(final Iterable<SAMRecord> noSecondaryAlignments,
                                                                                 final int gapSplitSensitivity,
                                                                                 final boolean splitGapped,
                                                                                 final Logger toolLogger) {

            Utils.validateArg(noSecondaryAlignments.iterator().hasNext(), "input collection of GATK reads is empty");

            final SAMRecord primaryAlignment
                    = Utils.stream(noSecondaryAlignments).filter(sam -> !sam.getSupplementaryAlignmentFlag())
                    .findFirst()
                    .orElseThrow(() -> new GATKException("no primary alignment for read " + noSecondaryAlignments.iterator().next().getReadName()));

            Utils.validate(!primaryAlignment.getCigar().containsOperator(CigarOperator.H),
                    "assumption that primary alignment does not contain hard clipping is invalid for read: " + primaryAlignment.toString());

            final byte[] contigSequence = primaryAlignment.getReadBases().clone();
            if (primaryAlignment.getReadNegativeStrandFlag()) {
                SequenceUtil.reverseComplement(contigSequence);
            }

            final Stream<AlignmentInterval> unSplitAIList = Utils.stream(noSecondaryAlignments).map(AlignmentInterval::new);
            final List<AlignmentInterval> parsedAlignments;
            if (splitGapped) {
                final int unClippedContigLength = primaryAlignment.getReadLength();
                parsedAlignments = unSplitAIList.map(ar ->
                        GappedAlignmentSplitter.split(ar, gapSplitSensitivity, unClippedContigLength))
                        .flatMap(Utils::stream).collect(Collectors.toList());
            } else {
                parsedAlignments = unSplitAIList.collect(Collectors.toList());
            }
            return new AlignedContig(primaryAlignment.getReadName(), contigSequence, parsedAlignments, false);
        }
    }

    /**
     * Makes sense out of the alignment records of the locally assembled contigs,
     * turn into annotated {@link VariantContext}'s, and writes them to VCF.
     */
    public static void discoverVariantsAndWriteVCF(final JavaRDD<AlignedContig> alignedContigs,
                                                   final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection parameters,
                                                   final Broadcast<ReferenceMultiSource> broadcastReference,
                                                   final String vcfOutputFileName,
                                                   final Logger localLogger,
                                                   final SAMSequenceDictionary samSequenceDictionary) {
        discoverVariantsAndWriteVCF(alignedContigs, parameters, broadcastReference, vcfOutputFileName,
                localLogger, null, null, samSequenceDictionary);
    }

    public static void discoverVariantsAndWriteVCF(final JavaRDD<AlignedContig> alignedContigs,
                                                   final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection parameters,
                                                   final Broadcast<ReferenceMultiSource> broadcastReference,
                                                   final String vcfOutputFileName, final Logger localLogger,
                                                   final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                   final ReadMetadata metadata,
                                                   final SAMSequenceDictionary sequenceDictionary) {

        Utils.validate(! (evidenceTargetLinks != null && metadata == null),
                "Must supply read metadata when incorporating evidence target links");

        JavaRDD<VariantContext> annotatedVariants =
                alignedContigs.filter(alignedContig -> alignedContig.alignmentIntervals.size()>1)                                     // filter out any contigs that has less than two alignment records
                        .mapToPair(alignedContig -> new Tuple2<>(alignedContig.contigSequence,                                        // filter a contig's alignment and massage into ordered collection of chimeric alignments
                                ChimericAlignment.parseOneContig(alignedContig, DEFAULT_MIN_ALIGNMENT_LENGTH, true, false)))
                        .flatMapToPair(DiscoverVariantsFromContigAlignmentsSAMSpark::discoverNovelAdjacencyFromChimericAlignments)    // a filter-passing contig's alignments may or may not produce novel adjacency
                        .groupByKey()                                                                                                 // group the same novel adjacency and alt haplotype sequence produced by different contigs together
                        .mapToPair(noveltyAndEvidence -> inferType(noveltyAndEvidence._1, noveltyAndEvidence._2))                     // type inference based on novel adjacency and evidence alignments
                        .map(noveltyTypeAndEvidence ->
                                annotateVariant(                                                                                      // annotate the novel adjacency and inferred type
                                        noveltyTypeAndEvidence._1.novelAdjacencyReferenceLocations,
                                        noveltyTypeAndEvidence._1.altHaplotypeSequence,
                                        noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2,
                                        broadcastReference));

        List<VariantContext> collectedAnnotatedVariants = annotatedVariants.collect();

        if (evidenceTargetLinks != null) {
            collectedAnnotatedVariants = processEvidenceTargetLinks(parameters,
                    localLogger,
                    evidenceTargetLinks,
                    metadata,
                    collectedAnnotatedVariants, broadcastReference.getValue());
        }

        SVVCFWriter.writeVCF(collectedAnnotatedVariants, vcfOutputFileName, sequenceDictionary, localLogger);
    }

    /**
     * Uses the input EvidenceTargetLinks to either annotate the variants called from assembly discovery with split
     * read and read pair evidence, or to call new imprecise variants if the number of pieces of evidence exceeds
     * a given threshold.
     */
    @VisibleForTesting
    static List<VariantContext> processEvidenceTargetLinks(final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection parameters,
                                                           final Logger localLogger,
                                                           final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                           final ReadMetadata metadata,
                                                           final List<VariantContext> assemblyDiscoveredVariants,
                                                           final ReferenceMultiSource reference) {
        final int originalEvidenceLinkSize = evidenceTargetLinks.size();
        List<VariantContext> result = assemblyDiscoveredVariants
                .stream()
                .map(variant -> AnnotatedVariantProducer.annotateWithImpreciseEvidenceLinks(
                        variant,
                        evidenceTargetLinks,
                        reference.getReferenceSequenceDictionary(null),
                        metadata, parameters.assemblyImpreciseEvidenceOverlapUncertainty))
                        .collect(Collectors.toList());
        localLogger.info("Used " + (originalEvidenceLinkSize - evidenceTargetLinks.size()) + " evidence target links to annotate assembled breakpoints");

        final List<VariantContext> impreciseVariants =
                Utils.stream(evidenceTargetLinks)
                        .map(p -> p._2)
                        .filter(EvidenceTargetLink::isImpreciseDeletion)
                        .filter(e -> e.getReadPairs() + e.getSplitReads() > parameters.impreciseEvidenceVariantCallingThreshold)
                        .map(e -> createImpreciseDeletionVariant(e, metadata, reference))
                        .collect(Collectors.toList());

        localLogger.info("Called " + impreciseVariants.size() + " imprecise deletion variants");
        result.addAll(impreciseVariants);
        return result;
    }

    private static VariantContext createImpreciseDeletionVariant(final EvidenceTargetLink e,
                                                                 final ReadMetadata metadata,
                                                                 final ReferenceMultiSource reference) {
        final SvType svType = new SimpleSVType.ImpreciseDeletion(e, metadata);
        return AnnotatedVariantProducer
                .produceAnnotatedVcFromEvidenceTargetLink(e, svType, metadata, reference);
    }

    // TODO: 7/6/17 interface to be changed in the new implementation, where one contig produces a set of NARL's.
    /**
     * Given contig alignments, emit novel adjacency not present on the reference to which the locally-assembled contigs were aligned.
     */
    public static Iterator<Tuple2<NovelAdjacencyAndInferredAltHaptype, ChimericAlignment>>
    discoverNovelAdjacencyFromChimericAlignments(final Tuple2<byte[], List<ChimericAlignment>> tigSeqAndChimeras) {
        return Utils.stream(tigSeqAndChimeras._2)
                .map(ca -> new Tuple2<>(new NovelAdjacencyReferenceLocations(ca, tigSeqAndChimeras._1), ca))
                .map(narlAndEvidenceCA -> {
                    final ChimericAlignment ca = narlAndEvidenceCA._2;
                    final SimpleInterval referenceSpan1 = ca.regionWithLowerCoordOnContig.referenceSpan;
                    final SimpleInterval referenceSpan2 = ca.regionWithHigherCoordOnContig.referenceSpan;
                    final byte[] altHaplotypeSeq;
                    if (ca.strandSwitch == StrandSwitch.NO_SWITCH && (referenceSpan1.contains(referenceSpan2) || referenceSpan2.contains(referenceSpan1))) {
                        altHaplotypeSeq = BreakpointComplications.extractAltHaplotypeForTandupExpansionWithContainment(ca.regionWithLowerCoordOnContig,
                                ca.regionWithHigherCoordOnContig, narlAndEvidenceCA._1.complication, tigSeqAndChimeras._1);
                    } else {
                        altHaplotypeSeq = null;
                    }
                    return new Tuple2<>(new NovelAdjacencyAndInferredAltHaptype(narlAndEvidenceCA._1, altHaplotypeSeq), ca);
                }).iterator();
    }

    // TODO: 7/6/17 interface to be changed in the new implementation, where a set of NRAL's associated with a single contig is considered together.
    /**
     * Given input novel adjacency and evidence chimeric alignments, infer type of variant.
     */
    public static Tuple2<NovelAdjacencyAndInferredAltHaptype, Tuple2<SvType, Iterable<ChimericAlignment>>> inferType(
            final NovelAdjacencyAndInferredAltHaptype novelAdjacency,
            final Iterable<ChimericAlignment> chimericAlignments) {
        return new Tuple2<>(novelAdjacency,
                new Tuple2<>(SvTypeInference.inferFromNovelAdjacency(novelAdjacency.novelAdjacencyReferenceLocations),
                        chimericAlignments));
    }

    /**
     * Produces annotated variant as described in {@link GATKSVVCFHeaderLines}, given input arguments.
     */
    public static VariantContext annotateVariant(final NovelAdjacencyReferenceLocations novelAdjacency,
                                                 final byte[] altHaplotypeSeq, final SvType inferredType,
                                                 final Iterable<ChimericAlignment> chimericAlignments,
                                                 final Broadcast<ReferenceMultiSource> broadcastReference)
            throws IOException {
        return AnnotatedVariantProducer
                .produceAnnotatedVcFromInferredTypeAndRefLocations(novelAdjacency.leftJustifiedLeftRefLoc,
                        novelAdjacency.leftJustifiedRightRefLoc.getStart(), novelAdjacency.complication,
                        inferredType, altHaplotypeSeq, chimericAlignments,
                        broadcastReference);
    }

    @DefaultSerializer(NovelAdjacencyAndInferredAltHaptype.Serializer.class)
    public static final class NovelAdjacencyAndInferredAltHaptype {

        private static final NovelAdjacencyReferenceLocations.Serializer localSerializer = new NovelAdjacencyReferenceLocations.Serializer();
        public final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations;
        public final byte[] altHaplotypeSequence;

        NovelAdjacencyAndInferredAltHaptype(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations,
                                            final byte[] altHaplotypeSequence) {
            this.novelAdjacencyReferenceLocations = novelAdjacencyReferenceLocations;
            this.altHaplotypeSequence = altHaplotypeSequence;
        }

        public NovelAdjacencyAndInferredAltHaptype(final Kryo kryo, final Input input) {
            novelAdjacencyReferenceLocations = localSerializer.read(kryo, input, NovelAdjacencyReferenceLocations.class);
            final boolean altSeqIsNull = input.readBoolean();
            if (altSeqIsNull) {
                altHaplotypeSequence = null;
            } else {
                final int arraySize = input.readInt();
                altHaplotypeSequence = new byte[arraySize];
                for (int i = 0 ; i < arraySize; ++i) {
                    altHaplotypeSequence[i] = input.readByte();
                }
            }
        }

        private void serialize(final Kryo kryo, final Output output) {
            localSerializer.write(kryo, output, novelAdjacencyReferenceLocations);
            if (altHaplotypeSequence==null) {
                output.writeBoolean(true);
            } else {
                output.writeBoolean(false);
                output.writeInt(altHaplotypeSequence.length);
                for (final byte b : altHaplotypeSequence) {
                    output.writeByte(b);
                }
            }
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<NovelAdjacencyAndInferredAltHaptype> {
            @Override
            public void write(final Kryo kryo, final Output output, final NovelAdjacencyAndInferredAltHaptype novelAdjacencyReferenceLocations ) {
                novelAdjacencyReferenceLocations.serialize(kryo, output);
            }

            @Override
            public NovelAdjacencyAndInferredAltHaptype read(final Kryo kryo, final Input input, final Class<NovelAdjacencyAndInferredAltHaptype> klass ) {
                return new NovelAdjacencyAndInferredAltHaptype(kryo, input);
            }
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            NovelAdjacencyAndInferredAltHaptype that = (NovelAdjacencyAndInferredAltHaptype) o;

            if (!novelAdjacencyReferenceLocations.equals(that.novelAdjacencyReferenceLocations)) return false;
            return Arrays.equals(altHaplotypeSequence, that.altHaplotypeSequence);
        }

        @Override
        public int hashCode() {
            int result = novelAdjacencyReferenceLocations.hashCode();
            result = 31 * result + Arrays.hashCode(altHaplotypeSequence);
            return result;
        }
    }
}
