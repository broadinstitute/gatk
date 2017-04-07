package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

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
        programGroup = SparkProgramGroup.class)
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

        final JavaRDD<AlignedAssembly.AlignedContig> parsedContigAlignments
                = new SAMFormattedContigAlignmentParser().loadAndFilterAndParseContigAlignments(ctx, localLogger, getReads(), getHeaderForReads());

        discoverVariantsAndWriteVCF(parsedContigAlignments, discoverStageArgs.fastaReference,
                ctx.broadcast(getReference()), getAuthenticatedGCSOptions(), vcfOutputFileName, localLogger);
    }

    public static final class SAMFormattedContigAlignmentParser implements AssemblyAlignmentParser {
        private static final long serialVersionUID = 1L;

        @SuppressWarnings("unchecked")
        public JavaRDD<AlignedAssembly.AlignedContig> loadAndFilterAndParseContigAlignments(final JavaSparkContext ctx,
                                                                                            final Logger toolLogger,
                                                                                            Object... localArguments) {
            final JavaRDD<GATKRead> unfilteredContigAlignments = (JavaRDD<GATKRead>) localArguments[0];
            final SAMFileHeader header = (SAMFileHeader) localArguments[1];

            return unfilteredContigAlignments
                    .filter(r -> !r.isSecondaryAlignment()).groupBy(GATKRead::getName).map(Tuple2::_2)
                    .map(iterable -> parseReadsAndBreakGaps(iterable, header, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, toolLogger));
        }

        /**
         * Iterates through the input {@code noSecondaryReads}, which are assumed to contain no secondary alignment (i.e. records with "XA" tag),
         * converts to custom {@link AlignedAssembly.AlignmentInterval} format and
         * break the records when the gap in the alignment reaches the specified {@code sensitivity}.
         * The size of the returned iterable of {@link AlignedAssembly.AlignmentInterval}'s is guaranteed to be no lower than that of the input iterable.
         */
        @VisibleForTesting
        public static AlignedAssembly.AlignedContig parseReadsAndBreakGaps(final Iterable<GATKRead> noSecondaryReads,
                                                                           final SAMFileHeader header,
                                                                           final int sensitivity,
                                                                           final Logger toolLogger) {

            Utils.validateArg(noSecondaryReads.iterator().hasNext(), "input collection of GATK reads is empty");

            final GATKRead primaryAlignment = Utils.stream(noSecondaryReads).filter(r -> !r.isSupplementaryAlignment())
                    .findFirst()
                    .orElseThrow(() -> new GATKException("no primary alignment for read " + noSecondaryReads.iterator().next().getName()));

            Utils.validate(!primaryAlignment.getCigar().containsOperator(CigarOperator.H),
                    "assumption that primary alignment does not contain hard clipping is invalid for read: " + primaryAlignment.toString());

            final byte[] contigSequence = primaryAlignment.getBases();
            if (primaryAlignment.isReverseStrand()) {
                SequenceUtil.reverseComplement(contigSequence);
            }
            final int unClippedContigLength = contigSequence.length;

            final List<AlignedAssembly.AlignmentInterval> parsedAlignments =
                    Utils.stream(noSecondaryReads)
                            .map(r -> AssemblyAlignmentParser.breakGappedAlignment(new AlignedAssembly.AlignmentInterval(r.convertToSAMRecord(header)), sensitivity, unClippedContigLength))
                            .flatMap(Utils::stream).collect(Collectors.toList());
            return new AlignedAssembly.AlignedContig(primaryAlignment.getName(), contigSequence, parsedAlignments);
        }
    }

    /**
     * Makes sense out of the alignment records of the locally assembled contigs,
     * turn into annotated {@link VariantContext}'s, and writes them to VCF.
     */
    public static void discoverVariantsAndWriteVCF(final JavaRDD<AlignedAssembly.AlignedContig> contigAlignments,
                                                   final String fastaReference, final Broadcast<ReferenceMultiSource> broadcastReference,
                                                   final PipelineOptions pipelineOptions, final String vcfFileName,
                                                   final Logger toolLogger) {

        final JavaRDD<VariantContext> variants
                = SVVariantConsensusDiscovery.discoverNovelAdjacencyFromChimericAlignments(contigAlignments, toolLogger)
                .map(tuple2 -> SVVariantConsensusDiscovery.discoverVariantsFromConsensus(tuple2, broadcastReference));

        SVVCFWriter.writeVCF(pipelineOptions, vcfFileName, fastaReference, variants, toolLogger);
    }
}
