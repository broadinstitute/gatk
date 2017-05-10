package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
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
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.util.List;
import java.util.stream.Collectors;

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
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs
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

        final JavaRDD<AlignedContig> parsedContigAlignments
                = new SAMFormattedContigAlignmentParser(getReads(), getHeaderForReads(), localLogger).getAlignedContigs();

        discoverVariantsAndWriteVCF(parsedContigAlignments, discoverStageArgs.fastaReference,
                ctx.broadcast(getReference()), getAuthenticatedGCSOptions(), vcfOutputFileName, localLogger);
    }

    public static final class SAMFormattedContigAlignmentParser extends AlignedContigGenerator implements Serializable {
        private static final long serialVersionUID = 1L;

        private final JavaRDD<GATKRead> unfilteredContigAlignments;
        private final SAMFileHeader header;
        private final Logger toolLogger;

        public SAMFormattedContigAlignmentParser(final JavaRDD<GATKRead> unfilteredContigAlignments,
                                                 final SAMFileHeader header,
                                                 final Logger toolLogger) {
            this.unfilteredContigAlignments = unfilteredContigAlignments;
            this.header = header;
            this.toolLogger = toolLogger;
        }

        @Override
        public JavaRDD<AlignedContig> getAlignedContigs() {
            return unfilteredContigAlignments
                    .filter(r -> !r.isSecondaryAlignment())
                    .groupBy(GATKRead::getName)
                    .map(Tuple2::_2)
                    .map(iterable -> parseReadsAndBreakGaps(Utils.stream(iterable).map(r->r.convertToSAMRecord(header)).collect(Collectors.toList()), header, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, toolLogger));
        }

        /**
         * Iterates through the input {@code noSecondaryAlignments}, which are assumed to contain no secondary alignment (i.e. records with "XA" tag),
         * converts to custom {@link AlignedAssembly.AlignmentInterval} format and
         * split the records when the gap in the alignment reaches the specified {@code sensitivity}.
         * The size of the returned iterable of {@link AlignedAssembly.AlignmentInterval}'s is guaranteed to be no lower than that of the input iterable.
         */
        @VisibleForTesting
        public static AlignedContig parseReadsAndBreakGaps(final Iterable<SAMRecord> noSecondaryAlignments,
                                                           final SAMFileHeader header,
                                                           final int sensitivity,
                                                           final Logger toolLogger) {

            Utils.validateArg(noSecondaryAlignments.iterator().hasNext(), "input collection of GATK reads is empty");

            final SAMRecord primaryAlignment = Utils.stream(noSecondaryAlignments).filter(sam -> !sam.getSupplementaryAlignmentFlag())
                    .findFirst()
                    .orElseThrow(() -> new GATKException("no primary alignment for read " + noSecondaryAlignments.iterator().next().getReadName()));

            Utils.validate(!primaryAlignment.getCigar().containsOperator(CigarOperator.H),
                    "assumption that primary alignment does not contain hard clipping is invalid for read: " + primaryAlignment.toString());

            final byte[] contigSequence = primaryAlignment.getReadBases().clone();
            if (primaryAlignment.getReadNegativeStrandFlag()) {
                SequenceUtil.reverseComplement(contigSequence);
            }
            final int unClippedContigLength = primaryAlignment.getReadLength();

            final List<AlignedAssembly.AlignmentInterval> parsedAlignments =
                    Utils.stream(noSecondaryAlignments)
                            .map(AlignedAssembly.AlignmentInterval::new)
                            .map(ar -> GappedAlignmentSplitter.split(ar, sensitivity, unClippedContigLength))
                            .flatMap(Utils::stream).collect(Collectors.toList());
            return new AlignedContig(primaryAlignment.getReadName(), contigSequence, parsedAlignments);
        }
    }

    /**
     * Makes sense out of the alignment records of the locally assembled contigs,
     * turn into annotated {@link VariantContext}'s, and writes them to VCF.
     */
    public static void discoverVariantsAndWriteVCF(final JavaRDD<AlignedContig> contigAlignments,
                                                   final String fastaReference, final Broadcast<ReferenceMultiSource> broadcastReference,
                                                   final PipelineOptions pipelineOptions, final String vcfFileName,
                                                   final Logger toolLogger) {

        final JavaRDD<VariantContext> variants
                = SVVariantConsensusDiscovery.discoverNovelAdjacencyFromChimericAlignments(contigAlignments, toolLogger)
                .map(tuple2 -> SVVariantConsensusDiscovery.discoverVariantsFromConsensus(tuple2, broadcastReference));

        SVVCFWriter.writeVCF(pipelineOptions, vcfFileName, fastaReference, variants, toolLogger);
    }
}
