package org.broadinstitute.hellbender.tools.spark.sv.discovery;

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
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.AlignmentSignatureBasicType.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.ReasonForAlignmentClassificationFailure.*;

/**
 * (Internal) Examines aligned contigs from local assemblies and calls structural variants or their breakpoints
 *
 * <p>
 *     This is an experimental tool and should not be of interest to most researchers. It is a prototype of a method
 *     for calling structural variation from alignments of assembled contigs and is under active development. For a
 *     more stable method for this, please see DiscoverVariantsFromContigAlignmentsSAMSpark.
 * </p>
 *
 * <p>
 *     This tool takes a file containing the alignments of assembled contigs (typically the output file produced by
 *     FindBreakpointEvidenceSpark) and searches for split alignments or alignments with large gaps indicating the
 *     presence of structural variation breakpoints. The type of each variation is determined by analyzing the
 *     signatures of the split alignments, and are written to VCF files in the designated output directory.
 * </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of assembled contigs aligned to reference.</li>
 *     <li>The reference to which the contigs have been aligned.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>Text files describing the discovered structural variants and complex structural variants in the specified output directory.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk SvDiscoverFromLocalAssemblyContigAlignmentsSpark \
 *     -I assemblies.sam \
 *     -R reference.2bit \
 *     -O output_directory
 * </pre>
 *
 * <h3>Notes</h3>
 * <p>The reference is broadcast by Spark, and must therefore be a .2bit file due to current restrictions.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Examines aligned contigs from local assemblies and calls structural variants or their breakpoints",
        summary =
        "This tool is used in development and should not be of interest to most researchers. It is a prototype of" +
        " structural variant calling, and has been under active developments. For more stable version," +
        " please see DiscoverVariantsFromContigAlignmentsSAMSpark." +
        " This tool takes a file containing the alignments of assembled contigs" +
        " (typically the output file produced by FindBreakpointEvidenceSpark) and searches for reads with" +
        " split alignments or large gaps indicating the presence of structural variation breakpoints." +
        " Variations' types are determined by analyzing the signatures of the split alignments," +
        " and are written to VCF files in the designated output directory.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class SvDiscoverFromLocalAssemblyContigAlignmentsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(SvDiscoverFromLocalAssemblyContigAlignmentsSpark.class);

    @ArgumentCollection
    private DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
            discoverStageArgs
            = new DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "file containing non-canonical chromosome names (e.g chrUn_KI270588v1) in the reference, human reference (hg19 or hg38) assumed when omitted",
            shortName = "alt-tigs",
            fullName = "non-canonical-contig-names-file", optional = true)
    public String nonCanonicalChromosomeNamesFile;

    @Argument(doc = "prefix for output files (including VCF files and if enabled, the signaling assembly contig's alignments); sample name will be appended after the provided argument",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputPrefix;

    @Argument(doc = "output SAM files", fullName = "write-sam", optional = true)
    private boolean writeSAMFiles;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.MAPPED);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast =
                StructuralVariationDiscoveryPipelineSpark.broadcastCNVCalls(ctx, getHeaderForReads(),
                        discoverStageArgs.cnvCallsFile);
        final String outputPrefixWithSampleName = outputPrefix + SVUtils.getSampleId(getHeaderForReads()) + "_";
        final SvDiscoveryInputMetaData svDiscoveryInputMetaData =
                new SvDiscoveryInputMetaData(ctx, discoverStageArgs, nonCanonicalChromosomeNamesFile, outputPrefixWithSampleName,
                        null, null, null,
                        cnvCallsBroadcast,
                        getHeaderForReads(), getReference(), localLogger);

        final AssemblyContigsClassifiedByAlignmentSignatures contigsByPossibleRawTypes =
                preprocess(svDiscoveryInputMetaData, getReads(), writeSAMFiles);

        dispatchJobs(contigsByPossibleRawTypes, svDiscoveryInputMetaData);
    }

    //==================================================================================================================

    public static final class AssemblyContigsClassifiedByAlignmentSignatures {
        final JavaRDD<AssemblyContigWithFineTunedAlignments> ambiguous;
        final JavaRDD<AssemblyContigWithFineTunedAlignments> incomplete;
        final JavaRDD<AssemblyContigWithFineTunedAlignments> misassembly;
        final JavaRDD<AssemblyContigWithFineTunedAlignments> simple;
        final JavaRDD<AssemblyContigWithFineTunedAlignments> complex;

        private AssemblyContigsClassifiedByAlignmentSignatures(final JavaRDD<AssemblyContigWithFineTunedAlignments> ambiguous,
                                                               final JavaRDD<AssemblyContigWithFineTunedAlignments> incomplete,
                                                               final JavaRDD<AssemblyContigWithFineTunedAlignments> misassembly,
                                                               final JavaRDD<AssemblyContigWithFineTunedAlignments> simple,
                                                               final JavaRDD<AssemblyContigWithFineTunedAlignments> complex) {
            this.ambiguous = ambiguous;
            this.incomplete = incomplete;
            this.misassembly = misassembly;
            this.simple = simple;
            this.complex = complex;
        }

        /**
         * Write SAM file, if requested, for original alignments of contigs recognized as "Ambiguous", "Incomplete", and "MisAssemblySuspect"
         * TODO: 11/17/17 salvation on assembly contigs that 1) has ambiguous "best" configuration, and 2) has incomplete picture; and flag accordingly
         *
         * TODO: 3/4/18 a bug is present here that even though only one alignment has not-bad MQ, it could contain a large gap,
         *      depending on the behavior of the other gap-less bad mappings,
         *      we may end up classifying the whole contig as incomplete, or signal-less,
         *      we should keep the single not-bad mapping and annotate accordingly
         */
        private void writeSAMfilesForSuspicious(final String outputPrefix, final JavaRDD<GATKRead> assemblyRawAlignments,
                                                final SAMFileHeader header) {
            final Set<String> ambiguousContigNames = new HashSet<>( ambiguous.map(AssemblyContigWithFineTunedAlignments::getContigName).collect() );
            final Set<String> incompleteContigNames = new HashSet<>( incomplete.map(AssemblyContigWithFineTunedAlignments::getContigName).collect() );
            final Set<String> misAssemblySuspectContigNames = new HashSet<>( misassembly.map(AssemblyContigWithFineTunedAlignments::getContigName).collect() );

            SvDiscoveryUtils.writeSAMRecords(assemblyRawAlignments, ambiguousContigNames, outputPrefix + AssemblyContigWithFineTunedAlignments.ReasonForAlignmentClassificationFailure.AMBIGUOUS.name() + ".bam", header);
            SvDiscoveryUtils.writeSAMRecords(assemblyRawAlignments, incompleteContigNames, outputPrefix + AssemblyContigWithFineTunedAlignments.ReasonForAlignmentClassificationFailure.INCOMPLETE.name() + ".bam", header);
            SvDiscoveryUtils.writeSAMRecords(assemblyRawAlignments, misAssemblySuspectContigNames, outputPrefix + AssemblyContigWithFineTunedAlignments.ReasonForAlignmentClassificationFailure.MIS_ASSEMBLY_OR_MAPPING_SUSPECT.name() + ".bam", header);
        }
    }

    /**
     * First parse the input alignments, then classify the assembly contigs based on their alignment signatures,
     * and return the contigs that are classified together for downstream inference.
     */
    public static AssemblyContigsClassifiedByAlignmentSignatures preprocess(final SvDiscoveryInputMetaData svDiscoveryInputMetaData,
                                                                            final JavaRDD<GATKRead> assemblyRawAlignments,
                                                                            final boolean writeSAMFiles) {

        final Broadcast<SAMFileHeader> headerBroadcast = svDiscoveryInputMetaData.sampleSpecificData.headerBroadcast;
        final Broadcast<Set<String>> canonicalChromosomesBroadcast = svDiscoveryInputMetaData.referenceData.canonicalChromosomesBroadcast;
        final Logger toolLogger = svDiscoveryInputMetaData.toolLogger;

        // filter alignments and split the gaps, hence the name "reconstructed"
        final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithChimericAlignmentsReconstructed =
                AssemblyContigAlignmentsConfigPicker
                        .createOptimalCoverageAlignmentSetsForContigs(assemblyRawAlignments, headerBroadcast.getValue(),
                                canonicalChromosomesBroadcast.getValue(), 0.0, toolLogger)
                        .filter(AssemblyContigWithFineTunedAlignments::isInformative).cache();
        toolLogger.info( contigsWithChimericAlignmentsReconstructed.count() +
                " contigs with chimeric alignments potentially giving SV signals.");

        // get assembly contigs by their possible type of SV based on studying alignment signature
        final AssemblyContigsClassifiedByAlignmentSignatures assemblyContigsClassifiedByAlignmentSignatures =
                classifyContigs(contigsWithChimericAlignmentsReconstructed, toolLogger);

        if (writeSAMFiles) {
            final String outputPrefix = svDiscoveryInputMetaData.outputPath;
            assemblyContigsClassifiedByAlignmentSignatures.writeSAMfilesForSuspicious(outputPrefix, assemblyRawAlignments,
                    headerBroadcast.getValue());
        }

        return assemblyContigsClassifiedByAlignmentSignatures;
    }

    private static AssemblyContigsClassifiedByAlignmentSignatures classifyContigs(final JavaRDD<AssemblyContigWithFineTunedAlignments> contigs,
                                                                                  final Logger toolLogger) {

        final JavaRDD<AssemblyContigWithFineTunedAlignments> suspiciousContigs = contigs
                .filter(tig ->
                        tig.getAlignmentSignatureBasicType()
                                .equals(UNKNOWN));
        final JavaRDD<AssemblyContigWithFineTunedAlignments> ambiguous = suspiciousContigs.filter(tig -> tig.getReasonForAlignmentClassificationFailure().equals(AMBIGUOUS));
        JavaRDD<AssemblyContigWithFineTunedAlignments> incomplete = suspiciousContigs.filter(tig -> tig.getReasonForAlignmentClassificationFailure().equals(INCOMPLETE));
        JavaRDD<AssemblyContigWithFineTunedAlignments> misassembly = suspiciousContigs.filter(tig -> tig.getReasonForAlignmentClassificationFailure().equals(MIS_ASSEMBLY_OR_MAPPING_SUSPECT));
        JavaRDD<AssemblyContigWithFineTunedAlignments> simple = contigs.filter(tig -> tig.getAlignmentSignatureBasicType().equals(SIMPLE));
        JavaRDD<AssemblyContigWithFineTunedAlignments> complex = contigs.filter(tig -> tig.getAlignmentSignatureBasicType().equals(COMPLEX));

        return new AssemblyContigsClassifiedByAlignmentSignatures(ambiguous, incomplete, misassembly, simple, complex);
    }

    //==================================================================================================================

    /**
     * Sends assembly contigs classified based on their alignment signature to
     * a corresponding breakpoint location inference unit.
     *
     * Two VCF files will be output: {@link #outputPrefix}"NonComplex.vcf" and {@link #outputPrefix}"Complex.vcf".
     *
     * Note that contigs with alignment signature classified as
     * {@link AssemblyContigWithFineTunedAlignments.AlignmentSignatureBasicType#UNKNOWN}
     * currently DO NOT generate any VCF yet.
     */
    public static void dispatchJobs(final AssemblyContigsClassifiedByAlignmentSignatures contigsByPossibleRawTypes,
                                    final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {

        final String outputPrefixWithSampleName = svDiscoveryInputMetaData.outputPath;

        // TODO: 1/10/18 bring back read annotation, see ticket 4228
        forNonComplexVariants(contigsByPossibleRawTypes.simple, svDiscoveryInputMetaData);

        final List<VariantContext> complexVariants =
                CpxVariantInterpreter.inferCpxVariant(contigsByPossibleRawTypes.complex, svDiscoveryInputMetaData);

        svDiscoveryInputMetaData.updateOutputPath(outputPrefixWithSampleName + "Complex.vcf");
        SVVCFWriter.writeVCF(complexVariants, svDiscoveryInputMetaData.outputPath,
                svDiscoveryInputMetaData.referenceData.referenceSequenceDictionaryBroadcast.getValue(),
                svDiscoveryInputMetaData.toolLogger);
    }

    private static void forNonComplexVariants(final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithSimpleChimera,
                                              final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputMetaData.referenceData.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputMetaData.referenceData.referenceSequenceDictionaryBroadcast;
        final String sampleId = svDiscoveryInputMetaData.sampleSpecificData.sampleId;
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputMetaData.sampleSpecificData.cnvCallsBroadcast;
        final String outputPrefixWithSampleName = svDiscoveryInputMetaData.outputPath;

        svDiscoveryInputMetaData.updateOutputPath(outputPrefixWithSampleName + "NonComplex.vcf");

        final List<VariantContext> annotatedSimpleVariants =
                new SimpleNovelAdjacencyInterpreter()
                        .inferTypeFromSingleContigSimpleChimera(contigsWithSimpleChimera, svDiscoveryInputMetaData)
                        .flatMap(pair ->
                            getVariantContextIterator(pair, sampleId, referenceBroadcast,
                                    referenceSequenceDictionaryBroadcast, cnvCallsBroadcast)
                        )
                        .collect();

        SVVCFWriter.writeVCF(annotatedSimpleVariants, svDiscoveryInputMetaData.outputPath,
                referenceSequenceDictionaryBroadcast.getValue(), svDiscoveryInputMetaData.toolLogger);
    }

    /**
     * This implementation is the 1st step going towards allowing re-interpretation,
     * below we simply take the inferred type and turn it to a VC,
     * future implementation may integrate other types of evidence and re-interpret if necessary
     */
    private static Iterator<VariantContext> getVariantContextIterator(final Tuple2<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>> pair,
                                                                      final String sampleId,
                                                                      final Broadcast<ReferenceMultiSource> referenceBroadcast,
                                                                      final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast,
                                                                      final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast)
            throws IOException {
        final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence = pair._1;
        final List<SvType> svTypes = pair._2;
        if( svTypes.isEmpty() || svTypes.size() > 2 ) {
            throw new GATKException("Wrong number of variants sent for analysis: " + pair._2.toString() +
                    "\nWe currently only support 1 (symbolic simple or CPX) or 2 (BND mate pairs) variants for producing annotated variants.");
        }
        if ( ! svTypes.get(0).isBreakEndOnly() ) { // simple SV type
            final NovelAdjacencyAndAltHaplotype narl = simpleNovelAdjacencyAndChimericAlignmentEvidence.getNovelAdjacencyReferenceLocations();
            final List<SimpleChimera> contigEvidence = simpleNovelAdjacencyAndChimericAlignmentEvidence.getAlignmentEvidence();

            if ( svTypes.size() == 2 ) { // RPL case with both path >= 50 bp
                final SvType firstVar = svTypes.get(0);
                final SvType secondVar = svTypes.get(1);
                final Tuple2<SvType, SvType> linkedVariants = new Tuple2<>(firstVar, secondVar);
                return AnnotatedVariantProducer.produceAnnotatedAndLinkedVcFromNovelAdjacency(linkedVariants,
                        simpleNovelAdjacencyAndChimericAlignmentEvidence,
                        referenceBroadcast, referenceSequenceDictionaryBroadcast, cnvCallsBroadcast, sampleId,
                        GATKSVVCFConstants.LINK).iterator();
            } else {
                final SvType inferredType = svTypes.get(0);

                final VariantContext variantContext = AnnotatedVariantProducer
                        .produceAnnotatedVcFromInferredTypeAndRefLocations(
                                narl, inferredType, contigEvidence,
                                referenceBroadcast, referenceSequenceDictionaryBroadcast, cnvCallsBroadcast, sampleId);
                return Collections.singletonList(variantContext).iterator();
            }
        } else { // BND mate pair
            final BreakEndVariantType firstMate = (BreakEndVariantType) svTypes.get(0);
            final BreakEndVariantType secondMate = (BreakEndVariantType) svTypes.get(1);

            final Tuple2<SvType, SvType> bndMates = new Tuple2<>(firstMate, secondMate);
            final List<VariantContext> variantContexts = AnnotatedVariantProducer
                    .produceAnnotatedAndLinkedVcFromNovelAdjacency(
                            bndMates, simpleNovelAdjacencyAndChimericAlignmentEvidence,
                            referenceBroadcast, referenceSequenceDictionaryBroadcast, cnvCallsBroadcast, sampleId,
                            GATKSVVCFConstants.BND_MATEID_STR);
            return variantContexts.iterator();
        }
    }

    //==================================================================================================================

    public static final class SAMFormattedContigAlignmentParser extends AlignedContigGenerator implements Serializable {
        private static final long serialVersionUID = 1L;

        private final JavaRDD<GATKRead> unfilteredContigAlignments;
        private final SAMFileHeader header;
        private final boolean splitGapped;

        public SAMFormattedContigAlignmentParser(final JavaRDD<GATKRead> unfilteredContigAlignments,
                                                 final SAMFileHeader header, final boolean splitGapped) {
            this.unfilteredContigAlignments = unfilteredContigAlignments;
            this.header = header;
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
                                    GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, splitGapped));
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
                                                                                 final boolean splitGapped) {

            Utils.validateArg(noSecondaryAlignments.iterator().hasNext(), "input collection of GATK reads is empty");

            final SAMRecord primaryAlignment
                    = Utils.stream(noSecondaryAlignments).filter(sam -> !sam.getSupplementaryAlignmentFlag())
                    .findFirst()
                    .orElseThrow(() -> new GATKException("no primary alignment for read " + noSecondaryAlignments.iterator().next().getReadName()));

            Utils.validate(!primaryAlignment.getCigar().containsOperator(CigarOperator.H),
                    "assumption that primary alignment does not contain hard clipping is invalid for read: " + primaryAlignment.toString());

            final byte[] contigSequence = primaryAlignment.getReadBases().clone();
            final List<AlignmentInterval> parsedAlignments;
            if ( primaryAlignment.getReadUnmappedFlag() ) { // the Cigar
                parsedAlignments = Collections.emptyList();
            } else {
                if (primaryAlignment.getReadNegativeStrandFlag()) {
                    SequenceUtil.reverseComplement(contigSequence);
                }

                final Stream<AlignmentInterval> unSplitAIList = Utils.stream(noSecondaryAlignments).map(AlignmentInterval::new);
                if (splitGapped) {
                    final int unClippedContigLength = primaryAlignment.getReadLength();
                    parsedAlignments = unSplitAIList.map(ar ->
                            ContigAlignmentsModifier.splitGappedAlignment(ar, gapSplitSensitivity, unClippedContigLength))
                            .flatMap(Utils::stream).collect(Collectors.toList());
                } else {
                    parsedAlignments = unSplitAIList.collect(Collectors.toList());
                }
            }
            return new AlignedContig(primaryAlignment.getReadName(), contigSequence, parsedAlignments);
        }
    }
}
