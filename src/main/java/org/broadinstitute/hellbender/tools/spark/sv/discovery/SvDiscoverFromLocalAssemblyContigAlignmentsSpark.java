package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
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
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.CpxVariantInterpreter;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleNovelAdjacencyInterpreter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFileUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.AlignmentSignatureBasicType.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.ReasonForAlignmentClassificationFailure;

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
    private String nonCanonicalChromosomeNamesFile;

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
        final String outputPrefixWithSampleName = getOutputPrefix();
        final SvDiscoveryInputMetaData svDiscoveryInputMetaData =
                new SvDiscoveryInputMetaData(ctx, discoverStageArgs, nonCanonicalChromosomeNamesFile, outputPrefixWithSampleName,
                        null, null, null,
                        cnvCallsBroadcast,
                        getHeaderForReads(), getReference(), localLogger);
        final JavaRDD<GATKRead> assemblyRawAlignments = getReads();

        final AssemblyContigsClassifiedByAlignmentSignatures contigsByPossibleRawTypes =
                preprocess(svDiscoveryInputMetaData, assemblyRawAlignments);

        dispatchJobs(ctx, contigsByPossibleRawTypes, svDiscoveryInputMetaData, assemblyRawAlignments, writeSAMFiles);
    }

    /**
     * @return prefix of outputs, with {@link #outputPrefix} decorated with sample name and trailing underscore
     */
    private String getOutputPrefix() {
        if ( Files.exists(Paths.get(outputPrefix)) ) {
            if (Files.isDirectory(Paths.get(outputPrefix))) // existing directory
                return outputPrefix + (outputPrefix.endsWith("/") ? "" : "/") + SVUtils.getSampleId(getHeaderForReads()) + "_";
            else
                throw new UserException("Provided prefix for output is pointing to an existing file: " + outputPrefix); // to avoid accidental override of a file
        } else { // prefixForOutput doesn't point to an existing file or directory
            return outputPrefix + (outputPrefix.endsWith("/") ? "" : "_") + SVUtils.getSampleId(getHeaderForReads()) + "_";
        }
    }

    //==================================================================================================================

    public static final class AssemblyContigsClassifiedByAlignmentSignatures {
        private final JavaRDD<AssemblyContigWithFineTunedAlignments> unknown;
        private final JavaRDD<AssemblyContigWithFineTunedAlignments> simple;
        private final JavaRDD<AssemblyContigWithFineTunedAlignments> complex;

        private AssemblyContigsClassifiedByAlignmentSignatures(final JavaRDD<AssemblyContigWithFineTunedAlignments> contigs) {
            unknown = contigs.filter(tig -> tig.getAlignmentSignatureBasicType().equals(UNKNOWN)).cache();
            simple = contigs.filter(tig -> tig.getAlignmentSignatureBasicType().equals(SIMPLE_CHIMERA)).cache();
            complex = contigs.filter(tig -> tig.getAlignmentSignatureBasicType().equals(COMPLEX)).cache();
        }

        public JavaRDD<AssemblyContigWithFineTunedAlignments> getContigsWithSignatureClassifiedAsUnknown() {
            return unknown;
        }

        public JavaRDD<AssemblyContigWithFineTunedAlignments> getContigsWithSignatureClassifiedAsSimpleChimera() {
            return simple;
        }

        public JavaRDD<AssemblyContigWithFineTunedAlignments> getContigsWithSignatureClassifiedAsComplex() {
            return complex;
        }

        /**
         * Write SAM file, if requested, for original alignments of contigs recognized as "Ambiguous", "Incomplete", and "MisAssemblySuspect"
         * TODO: 11/17/17 salvation on assembly contigs that 1) has ambiguous "best" configuration, and 2) has incomplete picture; and flag accordingly
         */
        private void writeSAMfilesForUnknown(final String outputPrefix, final JavaRDD<GATKRead> assemblyRawAlignments,
                                             final SAMFileHeader header) {

            final Map<String, ReasonForAlignmentClassificationFailure> tigNameToReason =
                    unknown.mapToPair(tig -> new Tuple2<>(tig.getContigName(), tig.getReasonForAlignmentClassificationFailure())).collectAsMap();

            final Set<String> namesOfInterest = new HashSet<>(tigNameToReason.keySet());

            final List<GATKRead> contigRawAlignments = assemblyRawAlignments
                    .filter(read -> namesOfInterest.contains(read.getName())).collect();

            final EnumMap<ReasonForAlignmentClassificationFailure, SAMFileWriter> writerForEachCase = new EnumMap<>(ReasonForAlignmentClassificationFailure.class);
            writerForEachCase.put(ReasonForAlignmentClassificationFailure.AMBIGUOUS,
                    SVFileUtils.getSAMFileWriter(outputPrefix + ReasonForAlignmentClassificationFailure.AMBIGUOUS.name() + ".bam",
                            header, false));
            writerForEachCase.put(ReasonForAlignmentClassificationFailure.INCOMPLETE,
                    SVFileUtils.getSAMFileWriter(outputPrefix + ReasonForAlignmentClassificationFailure.INCOMPLETE.name() + ".bam",
                            header, false));
            writerForEachCase.put(ReasonForAlignmentClassificationFailure.UNINFORMATIVE,
                    SVFileUtils.getSAMFileWriter(outputPrefix + ReasonForAlignmentClassificationFailure.UNINFORMATIVE.name() + ".bam",
                            header, false));
            contigRawAlignments.forEach(read -> {
                final ReasonForAlignmentClassificationFailure reason = tigNameToReason.get(read.getName());
                writerForEachCase.get(reason).addAlignment(read.convertToSAMRecord(header));
            });
        }
    }

    /**
     * First parse the input alignments, then classify the assembly contigs based on their alignment signatures,
     * and return the contigs that are classified together for downstream inference.
     */
    public static AssemblyContigsClassifiedByAlignmentSignatures preprocess(final SvDiscoveryInputMetaData svDiscoveryInputMetaData,
                                                                            final JavaRDD<GATKRead> assemblyRawAlignments) {

        final Broadcast<SAMFileHeader> headerBroadcast = svDiscoveryInputMetaData.getSampleSpecificData().getHeaderBroadcast();
        final Broadcast<Set<String>> canonicalChromosomesBroadcast = svDiscoveryInputMetaData.getReferenceData().getCanonicalChromosomesBroadcast();
        final Logger toolLogger = svDiscoveryInputMetaData.getToolLogger();

        final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithChimericAlignmentsReconstructed =
                AssemblyContigAlignmentsConfigPicker
                        .createOptimalCoverageAlignmentSetsForContigs(assemblyRawAlignments, headerBroadcast.getValue(),
                                canonicalChromosomesBroadcast.getValue(), 0.0, toolLogger)
                        .cache();
        toolLogger.info( contigsWithChimericAlignmentsReconstructed.count() +
                " contigs with chimeric alignments potentially giving SV signals.");

        return new AssemblyContigsClassifiedByAlignmentSignatures(contigsWithChimericAlignmentsReconstructed);
    }

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
    public static void dispatchJobs(final JavaSparkContext ctx,
                                    final AssemblyContigsClassifiedByAlignmentSignatures contigsByPossibleRawTypes,
                                    final SvDiscoveryInputMetaData svDiscoveryInputMetaData,
                                    final JavaRDD<GATKRead> assemblyRawAlignments,
                                    final boolean writeSAMFiles) {

        final String outputPrefixWithSampleName = svDiscoveryInputMetaData.getOutputPath();

        // TODO: 1/10/18 bring back read annotation, see ticket 4228

        final List<VariantContext> simpleVariants =
                SimpleNovelAdjacencyInterpreter.makeInterpretation(contigsByPossibleRawTypes.simple, svDiscoveryInputMetaData);
        contigsByPossibleRawTypes.simple.unpersist();
        SVVCFWriter.writeVCF(simpleVariants, outputPrefixWithSampleName + "NonComplex.vcf",
                svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast().getValue(),
                svDiscoveryInputMetaData.getToolLogger());

        final List<VariantContext> complexVariants =
                CpxVariantInterpreter.makeInterpretation(contigsByPossibleRawTypes.complex, svDiscoveryInputMetaData);
        contigsByPossibleRawTypes.complex.unpersist();
        SVVCFWriter.writeVCF(complexVariants, outputPrefixWithSampleName + "Complex.vcf",
                svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast().getValue(),
                svDiscoveryInputMetaData.getToolLogger());

        if (writeSAMFiles) {
            contigsByPossibleRawTypes.writeSAMfilesForUnknown(outputPrefixWithSampleName, assemblyRawAlignments,
                    svDiscoveryInputMetaData.getSampleSpecificData().getHeaderBroadcast().getValue());
        }

        final JavaRDD<VariantContext> complexVariantsRDD = ctx.parallelize(complexVariants);
        final SegmentedCpxVariantSimpleVariantExtractor.ExtractedSimpleVariants reInterpretedSimple =
                SegmentedCpxVariantSimpleVariantExtractor.extract(complexVariantsRDD, svDiscoveryInputMetaData, assemblyRawAlignments);
        final SAMSequenceDictionary refSeqDict = svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast().getValue();
        final Logger toolLogger = svDiscoveryInputMetaData.getToolLogger();
        final String derivedOneSegmentSimpleVCF = outputPrefixWithSampleName + "cpx_reinterpreted_simple_1_seg.vcf";
        final String derivedMultiSegmentSimpleVCF = outputPrefixWithSampleName + "cpx_reinterpreted_simple_multi_seg.vcf";
        SVVCFWriter.writeVCF(reInterpretedSimple.getReInterpretZeroOrOneSegmentCalls(), derivedOneSegmentSimpleVCF, refSeqDict, toolLogger);
        SVVCFWriter.writeVCF(reInterpretedSimple.getReInterpretMultiSegmentsCalls(), derivedMultiSegmentSimpleVCF, refSeqDict, toolLogger);
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
