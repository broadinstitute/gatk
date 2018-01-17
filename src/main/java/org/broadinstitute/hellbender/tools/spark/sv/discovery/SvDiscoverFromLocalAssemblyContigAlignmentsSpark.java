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
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFileUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
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
 * (Internal) Examines aligned contigs from local assemblies and calls complex structural variants
 *
 * <p>This tool is used in development and should not be of interest to most researchers.  It is a prototype of
 * complex structural variant calling, and has been superseded by other tools.</p>
 * <p>This tool takes a file containing the alignments of assembled contigs and searches for reads with
 * split alignments indicating the presence of structural variation breakpoints. The alignment signatures of the
 * split alignments are analyzed to determine the type of structural variation and written to a VCF file.</p>
 * <p>The input file is typically the output file produced by FindBreakpointEvidenceSpark.</p>
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
        oneLineSummary = "(Internal) Examines aligned contigs from local assemblies and calls complex structural variants",
        summary =
        "This tool is used in development and should not be of interest to most researchers.  It is a prototype of" +
        " complex structural variant calling, and has been superseded by other tools." +
        " This tool takes a file containing the alignments of assembled contigs and searches for reads with" +
        " split alignments indicating the presence of structural variation breakpoints. The alignment signatures of the" +
        " split alignments are analyzed to determine the type of structural variation and written to a VCF file." +
        " The input file is typically the output file produced by FindBreakpointEvidenceSpark.",
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

    @Argument(doc = "output directory for outputting VCF files for each type of variant, and if enabled, the signaling assembly contig's alignments",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDir;

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
    public ReadFilter makeReadFilter() {
        return ReadFilterLibrary.MAPPED;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast =
                StructuralVariationDiscoveryPipelineSpark.broadcastCNVCalls(ctx, getHeaderForReads(),
                        discoverStageArgs.cnvCallsFile);
        final SvDiscoveryInputData svDiscoveryInputData =
                new SvDiscoveryInputData(ctx, discoverStageArgs, outputDir,
                        null, null, null,
                        cnvCallsBroadcast,
                        getReads(), getHeaderForReads(), getReference(), localLogger);

        final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes =
                preprocess(svDiscoveryInputData, nonCanonicalChromosomeNamesFile, writeSAMFiles);

        dispatchJobs(contigsByPossibleRawTypes, svDiscoveryInputData);
    }

    //==================================================================================================================

    /**
     * First parse the input alignments, then classify the assembly contigs based on their alignment signatures,
     * and return the contigs that are classified together for downstream inference.
     */
    public static EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> preprocess(final SvDiscoveryInputData svDiscoveryInputData,
                                                                                               final String nonCanonicalChromosomeNamesFile,
                                                                                               final boolean writeSAMFiles) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final Broadcast<SAMFileHeader> headerBroadcast = svDiscoveryInputData.headerBroadcast;
        final JavaRDD<GATKRead> assemblyRawAlignments = svDiscoveryInputData.assemblyRawAlignments;
        final String outputPath = svDiscoveryInputData.outputPath;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        // filter alignments and split the gaps, hence the name "reconstructed"
        final JavaRDD<AlignedContig> contigsWithChimericAlignmentsReconstructed =
                AssemblyContigAlignmentsConfigPicker
                        .createOptimalCoverageAlignmentSetsForContigs(assemblyRawAlignments, headerBroadcast.getValue(),
                                nonCanonicalChromosomeNamesFile, 0.0, toolLogger)
                        .filter(lr -> lr.alignmentIntervals.size() > 1).cache();
        toolLogger.info( contigsWithChimericAlignmentsReconstructed.count() +
                " contigs with chimeric alignments potentially giving SV signals.");

        // classify assembly contigs by their possible type of SV based on studying alignment signature
        final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes =
                AssemblyContigAlignmentSignatureClassifier.classifyContigs(contigsWithChimericAlignmentsReconstructed,
                        referenceSequenceDictionaryBroadcast, toolLogger);

        try {
            IOUtils.createDirectory(outputPath);
            // write SAM file, if requested, for each type of possibly variant as recognized in {@link RawTypes}
            // and stored in {@code contigsByPossibleRawTypes} by extracting original alignments,
            if (writeSAMFiles) {
                contigsByPossibleRawTypes.forEach(
                        (k, v) ->
                                writeSAM(v, k.name(), assemblyRawAlignments, headerBroadcast, outputPath, toolLogger));
            }
        } catch (final IOException x) {
            throw new UserException.CouldNotCreateOutputFile("Could not create file at path:" +
                    outputPath + " due to " + x.getMessage(), x);
        }

        return contigsByPossibleRawTypes;
    }

    //==================================================================================================================

    /**
     * Sends assembly contigs classified based on their alignment signature to
     * a corresponding breakpoint location inference unit.
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#InsDel},
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#IntraChrStrandSwitch},
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#MappedInsertionBkpt}, and
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#Cpx} (PR to be reviewed and merged)
     * each will have its own VCF output in the directory specified in {@link #outputDir},
     * whereas
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#Incomplete},
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#Ambiguous}, and
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#MisAssemblySuspect}
     * currently DO NOT generate any VCF yet.
     * However, if flag {@link #writeSAMFiles} is turned on, alignments of all contigs that are classified to be any of
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes}
     * will be extracted and put in SAM files in {@link #outputDir} too.
     */
    public static void dispatchJobs(final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes,
                                    final SvDiscoveryInputData svDiscoveryInputData) {

        final String outputDir = svDiscoveryInputData.outputPath;

        // TODO: 1/10/18 bring back imprecise variants calling (definitely not here) and read annotation
        svDiscoveryInputData.updateOutputPath(outputDir+"/"+RawTypes.InsDel.name()+".vcf");
        new InsDelVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.InsDel), svDiscoveryInputData);

        svDiscoveryInputData.updateOutputPath(outputDir+"/"+RawTypes.IntraChrStrandSwitch.name()+".vcf");
        new SimpleStrandSwitchVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.IntraChrStrandSwitch), svDiscoveryInputData);

        svDiscoveryInputData.updateOutputPath(outputDir+"/"+RawTypes.MappedInsertionBkpt.name()+".vcf");
        new SuspectedTransLocDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.MappedInsertionBkpt), svDiscoveryInputData);

        svDiscoveryInputData.updateOutputPath(outputDir+"/"+RawTypes.Cpx.name()+".vcf");
        new CpxVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.Cpx), svDiscoveryInputData);
    }

    //==================================================================================================================

    /**
     * write SAM file for provided {@code filteredContigs}
     * by extracting original alignments from {@code originalAlignments},
     * to directory specified by {@code outputDir}.
     */
    private static void writeSAM(final JavaRDD<AssemblyContigWithFineTunedAlignments> filteredContigs, final String rawTypeString,
                                 final JavaRDD<GATKRead> originalAlignments, final Broadcast<SAMFileHeader> headerBroadcast,
                                 final String outputDir, final Logger toolLogger) {

        final Set<String> filteredReadNames = new HashSet<>( filteredContigs.map(decoratedTig -> decoratedTig.getSourceContig().contigName).distinct().collect() );
        toolLogger.info(filteredReadNames.size() + " contigs indicating " + rawTypeString);
        final JavaRDD<SAMRecord> splitLongReads = originalAlignments.filter(read -> filteredReadNames.contains(read.getName()))
                .map(read -> read.convertToSAMRecord(headerBroadcast.getValue()));
        SVFileUtils.writeSAMFile(outputDir+"/"+rawTypeString+".sam", splitLongReads.collect().iterator(),
                headerBroadcast.getValue(), false);
    }

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
            return new AlignedContig(primaryAlignment.getReadName(), contigSequence, parsedAlignments, false);
        }
    }
}
