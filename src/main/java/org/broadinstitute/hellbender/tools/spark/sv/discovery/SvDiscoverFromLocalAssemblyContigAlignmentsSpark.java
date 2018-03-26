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
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
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
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyContigAlignmentSignatureClassifier.RawTypes;

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
        final SvDiscoveryInputData svDiscoveryInputData =
                new SvDiscoveryInputData(ctx, discoverStageArgs, nonCanonicalChromosomeNamesFile, outputPrefixWithSampleName,
                        null, null, null,
                        cnvCallsBroadcast,
                        getReads(), getHeaderForReads(), getReference(), localLogger);

        final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes =
                preprocess(svDiscoveryInputData, writeSAMFiles);

        dispatchJobs(contigsByPossibleRawTypes, svDiscoveryInputData);
    }

    //==================================================================================================================

    /**
     * First parse the input alignments, then classify the assembly contigs based on their alignment signatures,
     * and return the contigs that are classified together for downstream inference.
     */
    public static EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> preprocess(final SvDiscoveryInputData svDiscoveryInputData,
                                                                                               final boolean writeSAMFiles) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceData.referenceSequenceDictionaryBroadcast;
        final Broadcast<SAMFileHeader> headerBroadcast = svDiscoveryInputData.sampleSpecificData.headerBroadcast;
        final Broadcast<Set<String>> canonicalChromosomesBroadcast = svDiscoveryInputData.referenceData.canonicalChromosomesBroadcast;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        final JavaRDD<GATKRead> assemblyRawAlignments = svDiscoveryInputData.sampleSpecificData.assemblyRawAlignments;

        // filter alignments and split the gaps, hence the name "reconstructed"
        final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithChimericAlignmentsReconstructed =
                AssemblyContigAlignmentsConfigPicker
                        .createOptimalCoverageAlignmentSetsForContigs(assemblyRawAlignments, headerBroadcast.getValue(),
                                canonicalChromosomesBroadcast.getValue(), 0.0, toolLogger)
                        .filter(AssemblyContigWithFineTunedAlignments::isInformative).cache();
        toolLogger.info( contigsWithChimericAlignmentsReconstructed.count() +
                " contigs with chimeric alignments potentially giving SV signals.");

        // classify assembly contigs by their possible type of SV based on studying alignment signature
        final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes =
                AssemblyContigAlignmentSignatureClassifier.classifyContigs(contigsWithChimericAlignmentsReconstructed,
                        referenceSequenceDictionaryBroadcast, toolLogger);

        // write SAM file, if requested, for original alignments of contigs recognized as "Ambiguous", "Incomplete", and "MisAssemblySuspect"
        if (writeSAMFiles) {
            final String outputPrefix = svDiscoveryInputData.outputPath;

            final Set<String> ambiguousContigNames = new HashSet<>( contigsByPossibleRawTypes.get(RawTypes.Ambiguous).map(AssemblyContigWithFineTunedAlignments::getContigName).collect() );
            final Set<String> incompleteContigNames = new HashSet<>( contigsByPossibleRawTypes.get(RawTypes.Incomplete).map(AssemblyContigWithFineTunedAlignments::getContigName).collect() );
            final Set<String> misAssemblySuspectContigNames = new HashSet<>( contigsByPossibleRawTypes.get(RawTypes.MisAssemblySuspect).map(AssemblyContigWithFineTunedAlignments::getContigName).collect() );

            final SAMFileHeader header = headerBroadcast.getValue();
            SvDiscoveryUtils.writeSAMRecords(assemblyRawAlignments, ambiguousContigNames, outputPrefix + RawTypes.Ambiguous.name() + ".bam", header);
            SvDiscoveryUtils.writeSAMRecords(assemblyRawAlignments, incompleteContigNames, outputPrefix + RawTypes.Incomplete.name() + ".bam", header);
            SvDiscoveryUtils.writeSAMRecords(assemblyRawAlignments, misAssemblySuspectContigNames, outputPrefix + RawTypes.MisAssemblySuspect.name() + ".bam", header);
        }

        return contigsByPossibleRawTypes;
    }

    //==================================================================================================================

    /**
     * Sends assembly contigs classified based on their alignment signature to
     * a corresponding breakpoint location inference unit.
     *
     * Two VCF files will be output: {@link #outputPrefix}"NonComplex.vcf" and {@link #outputPrefix}"Cpx.vcf".
     *
     * Note that
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#Incomplete},
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#Ambiguous}, and
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#MisAssemblySuspect}
     * currently DO NOT generate any VCF yet.
     * However, if flag {@link #writeSAMFiles} is turned on,
     * alignments of all contigs that are classified to be any of
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes}
     * will be extracted and put in SAM files in {@link #outputPrefix} too.
     */
    public static void dispatchJobs(final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes,
                                    final SvDiscoveryInputData svDiscoveryInputData) {

        final String outputPrefixWithSampleName = svDiscoveryInputData.outputPath;

        // TODO: 1/10/18 bring back read annotation, see ticket 4228
        forNonComplexVariants(contigsByPossibleRawTypes, svDiscoveryInputData);

        final List<VariantContext> complexVariants =
                CpxVariantInterpreter.inferCpxVariant(contigsByPossibleRawTypes.get(RawTypes.Cpx), svDiscoveryInputData);

        svDiscoveryInputData.updateOutputPath(outputPrefixWithSampleName + RawTypes.Cpx.name() + ".vcf");
        SVVCFWriter.writeVCF(complexVariants, svDiscoveryInputData.outputPath,
                svDiscoveryInputData.referenceData.referenceSequenceDictionaryBroadcast.getValue(),
                svDiscoveryInputData.toolLogger);
    }

    private static void forNonComplexVariants(final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes,
                                              final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceData.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceData.referenceSequenceDictionaryBroadcast;
        final String sampleId = svDiscoveryInputData.sampleSpecificData.sampleId;
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputData.sampleSpecificData.cnvCallsBroadcast;
        final String outputPrefixWithSampleName = svDiscoveryInputData.outputPath;

        svDiscoveryInputData.updateOutputPath(outputPrefixWithSampleName + "NonComplex.vcf");

        final JavaRDD<AssemblyContigWithFineTunedAlignments> nonComplexSignatures =
                contigsByPossibleRawTypes.get(RawTypes.InsDel)
                        .union(contigsByPossibleRawTypes.get(RawTypes.IntraChrStrandSwitch))
                        .union(contigsByPossibleRawTypes.get(RawTypes.MappedInsertionBkpt));

        final List<VariantContext> annotatedSimpleVariants =
                new SimpleNovelAdjacencyInterpreter()
                        .inferTypeFromSingleContigSimpleChimera(nonComplexSignatures, svDiscoveryInputData)
                        .flatMap(pair ->
                            getVariantContextIterator(pair, sampleId, referenceBroadcast,
                                    referenceSequenceDictionaryBroadcast, cnvCallsBroadcast)
                        )
                        .collect();

        SVVCFWriter.writeVCF(annotatedSimpleVariants, svDiscoveryInputData.outputPath,
                referenceSequenceDictionaryBroadcast.getValue(), svDiscoveryInputData.toolLogger);
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
        if (svTypes.size() == 1) { // simple SV type
            final SvType inferredType = svTypes.get(0);
            final NovelAdjacencyAndAltHaplotype narl = simpleNovelAdjacencyAndChimericAlignmentEvidence.getNovelAdjacencyReferenceLocations();
            final SimpleInterval variantPos = narl.getLeftJustifiedLeftRefLoc();
            final int end = narl.getLeftJustifiedRightRefLoc().getEnd();
            final VariantContext variantContext = AnnotatedVariantProducer
                    .produceAnnotatedVcFromInferredTypeAndRefLocations(variantPos, end, narl.getComplication(),
                            inferredType, simpleNovelAdjacencyAndChimericAlignmentEvidence.getAltHaplotypeSequence(),
                            simpleNovelAdjacencyAndChimericAlignmentEvidence.getAlignmentEvidence(),
                            referenceBroadcast, referenceSequenceDictionaryBroadcast, cnvCallsBroadcast, sampleId);
            return Collections.singletonList(variantContext).iterator();
        } else { // BND mate pair
            final BreakEndVariantType firstMate = (BreakEndVariantType) svTypes.get(0);
            final BreakEndVariantType secondMate = (BreakEndVariantType) svTypes.get(1);

            final Tuple2<BreakEndVariantType, BreakEndVariantType> bndMates = new Tuple2<>(firstMate, secondMate);
            final List<VariantContext> variantContexts = AnnotatedVariantProducer
                    .produceAnnotatedBNDmatesVcFromNovelAdjacency(
                            simpleNovelAdjacencyAndChimericAlignmentEvidence.getNovelAdjacencyReferenceLocations(),
                            bndMates,
                            simpleNovelAdjacencyAndChimericAlignmentEvidence.getAlignmentEvidence(),
                            referenceBroadcast, referenceSequenceDictionaryBroadcast, sampleId);
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
