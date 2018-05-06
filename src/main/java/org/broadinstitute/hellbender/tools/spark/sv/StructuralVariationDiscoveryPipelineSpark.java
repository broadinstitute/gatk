package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputMetaData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContigGenerator;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ImpreciseVariantDetector;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.FindBreakpointEvidenceSpark;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Serializable;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Runs the structural variation discovery workflow on a single sample
 *
 * <p>This tool packages the algorithms described in {@link FindBreakpointEvidenceSpark} and
 * {@link org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark}
 * as an integrated workflow.  Please consult the
 * descriptions of those tools for more details about the algorithms employed.  In brief, input reads are examined
 * for evidence of structural variation in a genomic region, regions so identified are locally assembled, and
 * the local assemblies are called for structural variation.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of aligned reads.</li>
 *     <li>The reference to which the reads have been aligned.</li>
 *     <li>A BWA index image for the reference.
 *         You can use BwaMemIndexImageCreator to create the index image file.</li>
 *     <li>A list of ubiquitous kmers to ignore.
 *         You can use FindBadGenomicGenomicKmersSpark to create the list of kmers to ignore.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A vcf file describing the discovered structural variants.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk StructuralVariationDiscoveryPipelineSpark \
 *     -I input_reads.bam \
 *     -R reference.2bit \
 *     --aligner-index-image reference.img \
 *     --kmers-to-ignore ignored_kmers.txt \
 *     -O structural_variants.vcf
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 *
 * <h3>Caveats</h3>
 * <p>Expected input is a paired-end, coordinate-sorted BAM with around 30x coverage.
 * Coverage much lower than that probably won't work well.</p>
 * <p>The reference is broadcast by Spark, and must therefore be a .2bit file due to current restrictions.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Runs the structural variation discovery workflow on a single sample",
        summary =
        "This tool packages the algorithms described in FindBreakpointEvidenceSpark and" +
        " DiscoverVariantsFromContigAlignmentsSAMSpark as an integrated workflow.  Please consult the" +
        " descriptions of those tools for more details about the algorithms employed.  In brief, input reads are examined" +
        " for evidence of structural variation in a genomic region, regions so identified are locally assembled, and" +
        " the local assemblies are called for structural variation.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class StructuralVariationDiscoveryPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(StructuralVariationDiscoveryPipelineSpark.class);
    @ArgumentCollection
    private final FindBreakpointEvidenceSparkArgumentCollection evidenceAndAssemblyArgs
            = new FindBreakpointEvidenceSparkArgumentCollection();
    @ArgumentCollection
    private final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs
            = new DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "sam file for aligned contigs", fullName = "contig-sam-file")
    private String outputAssemblyAlignments;

    @Argument(doc = "directory for VCF output, including those from experimental interpretation tool if so requested, " +
            "will be created if not present; sample name will be appended after the provided argument",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String variantsOutDir;

    @Advanced
    @Argument(doc = "flag to signal that user wants to run experimental interpretation tool as well",
            fullName = "exp-interpret", optional = true)
    private Boolean expInterpret = false;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    public boolean requiresReference() {return true;}

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        Utils.validate(evidenceAndAssemblyArgs.externalEvidenceFile == null || discoverStageArgs.cnvCallsFile == null,
                "Please only specify one of externalEvidenceFile or cnvCallsFile");

        if (discoverStageArgs.cnvCallsFile != null) {
            evidenceAndAssemblyArgs.externalEvidenceFile = discoverStageArgs.cnvCallsFile;
        }

        JavaRDD<GATKRead> unfilteredReads = getUnfilteredReads();
        final SAMFileHeader headerForReads = getHeaderForReads();

        // gather evidence, run assembly, and align
        final FindBreakpointEvidenceSpark.AssembledEvidenceResults assembledEvidenceResults =
                FindBreakpointEvidenceSpark
                        .gatherEvidenceAndWriteContigSamFile(ctx,
                                evidenceAndAssemblyArgs,
                                headerForReads,
                                unfilteredReads,
                                outputAssemblyAlignments,
                                localLogger);

        // todo: when we call imprecise variants don't return here
        if (assembledEvidenceResults.getAlignedAssemblyOrExcuseList().isEmpty()) return;

        // parse the contig alignments and extract necessary information
        final JavaRDD<AlignedContig> parsedAlignments =
                new InMemoryAlignmentParser(ctx, assembledEvidenceResults.getAlignedAssemblyOrExcuseList(), headerForReads)
                        .getAlignedContigs();

        // todo: when we call imprecise variants don't return here
        if(parsedAlignments.isEmpty()) return;

        final SvDiscoveryInputMetaData svDiscoveryInputMetaData = getSvDiscoveryInputData(ctx, headerForReads, assembledEvidenceResults);

        // TODO: 1/14/18 this is to be phased-out: old way of calling precise variants
        // assembled breakpoints
        @SuppressWarnings("deprecation")
        List<VariantContext> assemblyBasedVariants =
                org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark
                        .discoverVariantsFromChimeras(svDiscoveryInputMetaData, parsedAlignments);

        final List<VariantContext> annotatedVariants = processEvidenceTargetLinks(assemblyBasedVariants, svDiscoveryInputMetaData);

        final String outputPath = svDiscoveryInputMetaData.outputPath;
        final SAMSequenceDictionary refSeqDictionary = svDiscoveryInputMetaData.referenceData.referenceSequenceDictionaryBroadcast.getValue();
        final Logger toolLogger = svDiscoveryInputMetaData.toolLogger;
        SVVCFWriter.writeVCF(annotatedVariants, outputPath + "inv_del_ins.vcf", refSeqDictionary, toolLogger);

        // TODO: 1/14/18 this is the next version of precise variant calling
        if ( expInterpret != null ) {
            experimentalInterpretation(ctx, expInterpret, assembledEvidenceResults, svDiscoveryInputMetaData);
        }
    }

    private SvDiscoveryInputMetaData getSvDiscoveryInputData(final JavaSparkContext ctx,
                                                             final SAMFileHeader headerForReads,
                                                             final FindBreakpointEvidenceSpark.AssembledEvidenceResults assembledEvidenceResults) {
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast =
                broadcastCNVCalls(ctx, headerForReads, discoverStageArgs.cnvCallsFile);
        try {
            if ( !java.nio.file.Files.exists(Paths.get(variantsOutDir)) ) {
                IOUtils.createDirectory(variantsOutDir);
            }
        } catch (final IOException ioex) {
            throw new GATKException("Failed to create output directory " + variantsOutDir + " though it does not yet exist", ioex);
        }

        final String outputPrefixWithSampleName = variantsOutDir + (variantsOutDir.endsWith("/") ? "" : "/")
                                                    + SVUtils.getSampleId(headerForReads) + "_";

        return new SvDiscoveryInputMetaData(ctx, discoverStageArgs, evidenceAndAssemblyArgs.crossContigsToIgnoreFile,
                outputPrefixWithSampleName,
                assembledEvidenceResults.getReadMetadata(), assembledEvidenceResults.getAssembledIntervals(),
                makeEvidenceLinkTree(assembledEvidenceResults.getEvidenceTargetLinks()),
                cnvCallsBroadcast, getHeaderForReads(), getReference(), localLogger);
    }

    /**
     * Uses the input EvidenceTargetLinks to
     *  <ul>
     *      <li>
     *          either annotate the variants called from assembly discovery with split read and read pair evidence, or
     *      </li>
     *      <li>
     *          to call new imprecise variants if the number of pieces of evidence exceeds a given threshold.
     *      </li>
     *  </ul>
     *
     */
    private static List<VariantContext> processEvidenceTargetLinks(List<VariantContext> assemblyBasedVariants,
                                                                   final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {

        final List<VariantContext> annotatedVariants;
        if (svDiscoveryInputMetaData.sampleSpecificData.evidenceTargetLinks != null) {
            final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks = svDiscoveryInputMetaData.sampleSpecificData.evidenceTargetLinks;
            final ReadMetadata readMetadata = svDiscoveryInputMetaData.sampleSpecificData.readMetadata;
            final ReferenceMultiSource reference = svDiscoveryInputMetaData.referenceData.referenceBroadcast.getValue();
            final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputMetaData.discoverStageArgs;
            final Logger toolLogger = svDiscoveryInputMetaData.toolLogger;

            // annotate with evidence links
            annotatedVariants = AnnotatedVariantProducer.
                    annotateBreakpointBasedCallsWithImpreciseEvidenceLinks(assemblyBasedVariants,
                            evidenceTargetLinks, readMetadata, reference, discoverStageArgs, toolLogger);

            // then also imprecise deletion
            final List<VariantContext> impreciseVariants = ImpreciseVariantDetector.
                    callImpreciseDeletionFromEvidenceLinks(evidenceTargetLinks, readMetadata, reference,
                            discoverStageArgs.impreciseVariantEvidenceThreshold,
                            discoverStageArgs.maxCallableImpreciseVariantDeletionSize,
                            toolLogger);

            annotatedVariants.addAll(impreciseVariants);
        } else {
            annotatedVariants = assemblyBasedVariants;
        }

        return annotatedVariants;
    }

    private static void experimentalInterpretation(final JavaSparkContext ctx, final Boolean expInterpret,
                                                   final FindBreakpointEvidenceSpark.AssembledEvidenceResults assembledEvidenceResults,
                                                   final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {

        if ( !expInterpret)
            return;

        final JavaRDD<GATKRead> assemblyRawAlignments = getContigRawAlignments(ctx, assembledEvidenceResults, svDiscoveryInputMetaData);

        final String updatedOutputPath = svDiscoveryInputMetaData.outputPath + "experimentalInterpretation_";
        svDiscoveryInputMetaData.updateOutputPath(updatedOutputPath);

        SvDiscoverFromLocalAssemblyContigAlignmentsSpark.AssemblyContigsClassifiedByAlignmentSignatures contigsByPossibleRawTypes
                = SvDiscoverFromLocalAssemblyContigAlignmentsSpark.preprocess(svDiscoveryInputMetaData, assemblyRawAlignments);

        SvDiscoverFromLocalAssemblyContigAlignmentsSpark.dispatchJobs(contigsByPossibleRawTypes, svDiscoveryInputMetaData,
                assemblyRawAlignments, true);
    }

    private static JavaRDD<GATKRead> getContigRawAlignments(final JavaSparkContext ctx,
                                                            final FindBreakpointEvidenceSpark.AssembledEvidenceResults assembledEvidenceResults,
                                                            final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast =
                svDiscoveryInputMetaData.referenceData.referenceSequenceDictionaryBroadcast;
        final Broadcast<SAMFileHeader> headerBroadcast = svDiscoveryInputMetaData.sampleSpecificData.headerBroadcast;
        final SAMFileHeader headerForReads = headerBroadcast.getValue();
        final SAMReadGroupRecord contigAlignmentsReadGroup = new SAMReadGroupRecord(SVUtils.GATKSV_CONTIG_ALIGNMENTS_READ_GROUP_ID);
        final List<String> refNames = SequenceDictionaryUtils.getContigNamesList(referenceSequenceDictionaryBroadcast.getValue());

        return ctx.parallelize(
                assembledEvidenceResults
                        .getAlignedAssemblyOrExcuseList().stream()
                        .filter(AlignedAssemblyOrExcuse::isNotFailure)
                        .flatMap(aa -> aa.toSAMStreamForAlignmentsOfThisAssembly(headerForReads, refNames, contigAlignmentsReadGroup))
                        .map(SAMRecordToGATKReadAdapter::new)
                        .collect(Collectors.toList())
        );
    }

    /**
     * Makes a PairedStrandedIntervalTree from a list of EvidenceTargetLinks. The value of each entry in the resulting tree
     * will be the original EvidenceTargetLink. If the input list is null, returns a null tree.
     */
    private PairedStrandedIntervalTree<EvidenceTargetLink> makeEvidenceLinkTree(final List<EvidenceTargetLink> evidenceTargetLinks) {
        final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceLinkTree;

        if (evidenceTargetLinks != null) {
            evidenceLinkTree = new PairedStrandedIntervalTree<>();
            evidenceTargetLinks.forEach(l -> evidenceLinkTree.put(l.getPairedStrandedIntervals(), l));
        } else {
            evidenceLinkTree = null;
        }
        return evidenceLinkTree;
    }

    public static final class InMemoryAlignmentParser extends AlignedContigGenerator implements Serializable {
        private static final long serialVersionUID = 1L;

        private final JavaSparkContext ctx;
        private final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList;
        private final SAMFileHeader header;


        InMemoryAlignmentParser(final JavaSparkContext ctx,
                                final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList,
                                final SAMFileHeader header) {
            this.ctx = ctx;
            this.alignedAssemblyOrExcuseList = alignedAssemblyOrExcuseList;
            this.header = header;
        }

        /**
         * Filters out "failed" assemblies, unmapped and secondary (i.e. "XA") alignments, and
         * turn the alignments of contigs into custom {@link AlignmentInterval} format.
         * Should be generating the same output as {@link #filterAndConvertToAlignedContigDirect(Iterable, List, SAMFileHeader)};
         * and currently this assertion is tested in {@see AlignedContigGeneratorUnitTest#testConvertAlignedAssemblyOrExcuseToAlignedContigsDirectAndConcordanceWithSAMRoute()}
         */
        @VisibleForTesting
        public static JavaRDD<AlignedContig> filterAndConvertToAlignedContigViaSAM(final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList,
                                                                                   final SAMFileHeader header,
                                                                                   final JavaSparkContext ctx) {

            final SAMFileHeader cleanHeader = new SAMFileHeader(header.getSequenceDictionary());
            final List<String> refNames = SequenceDictionaryUtils.getContigNamesList(header.getSequenceDictionary());

            return ctx.parallelize(alignedAssemblyOrExcuseList)
                    .filter(AlignedAssemblyOrExcuse::isNotFailure)
                    .flatMap(alignedAssemblyNoExcuse -> {
                                final FermiLiteAssembly assembly = alignedAssemblyNoExcuse.getAssembly();
                                final int assemblyId = alignedAssemblyNoExcuse.getAssemblyId();
                                final List<List<BwaMemAlignment>> allAlignmentsOfThisAssembly = alignedAssemblyNoExcuse.getContigAlignments();
                                final int nContigs = assembly.getNContigs();
                                return IntStream.range(0, nContigs)
                                        .mapToObj(contigIdx ->
                                                BwaMemAlignmentUtils.toSAMStreamForRead(
                                                        AlignedAssemblyOrExcuse.formatContigName(assemblyId, contigIdx),
                                                        assembly.getContig(contigIdx).getSequence(),
                                                        allAlignmentsOfThisAssembly.get(contigIdx),
                                                        cleanHeader, refNames,
                                                        new SAMReadGroupRecord(SVUtils.GATKSV_CONTIG_ALIGNMENTS_READ_GROUP_ID)
                                                )
                                        ).iterator();
                            }
                    )
                    .map(forOneContig ->
                            forOneContig.filter(sam -> !sam.getReadUnmappedFlag() && !sam.isSecondaryAlignment())
                                    .collect(Collectors.toList()))
                    .filter(list -> !list.isEmpty())
                    .map(forOneContig ->
                            SvDiscoverFromLocalAssemblyContigAlignmentsSpark.
                                    SAMFormattedContigAlignmentParser.
                                    parseReadsAndOptionallySplitGappedAlignments(forOneContig,
                                            StructuralVariationDiscoveryArgumentCollection
                                                    .DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
                                                    .GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY,
                                            true));
        }

        /**
         * Filter out "failed" assemblies, unmapped and secondary (i.e. "XA") alignments, and
         * turn the alignments of contigs into custom {@link AlignmentInterval} format.
         * Should be generating the same output as {@link #filterAndConvertToAlignedContigViaSAM(List, SAMFileHeader, JavaSparkContext)};
         * and currently this assertion is tested in {@see AlignedContigGeneratorUnitTest#testConvertAlignedAssemblyOrExcuseToAlignedContigsDirectAndConcordanceWithSAMRoute()}
         */
        @VisibleForTesting
        public static List<AlignedContig> filterAndConvertToAlignedContigDirect(final Iterable<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseIterable,
                                                                                final List<String> refNames, final SAMFileHeader header) {

            return Utils.stream(alignedAssemblyOrExcuseIterable)
                    .filter(AlignedAssemblyOrExcuse::isNotFailure)
                    .map(alignedAssembly -> getAlignedContigsInOneAssembly(alignedAssembly, refNames, header))
                    .flatMap(Utils::stream)                                     // size == total # of contigs' from all successful assemblies
                    .filter(contig -> !contig.getAlignments().isEmpty())     // filter out unmapped and contigs without primary alignments
                    .collect(Collectors.toList());
        }

        /**
         * Work on "successful" assembly and turn its contigs' alignments to custom {@link AlignmentInterval} format.
         */
        @VisibleForTesting
        public static Iterable<AlignedContig> getAlignedContigsInOneAssembly(final AlignedAssemblyOrExcuse alignedAssembly,
                                                                             final List<String> refNames,
                                                                             final SAMFileHeader header) {

            final FermiLiteAssembly assembly = alignedAssembly.getAssembly();

            final List<List<BwaMemAlignment>> allAlignments = alignedAssembly.getContigAlignments();

            return IntStream.range(0, assembly.getNContigs())
                    .mapToObj( contigIdx -> {
                        final byte[] contigSequence = assembly.getContig(contigIdx).getSequence();
                        final String contigName = AlignedAssemblyOrExcuse.formatContigName(alignedAssembly.getAssemblyId(), contigIdx);
                        final List<AlignmentInterval> alignmentsForOneContig
                                = getAlignmentsForOneContig(contigName, contigSequence, allAlignments.get(contigIdx), refNames, header);
                        return new AlignedContig(contigName, contigSequence, alignmentsForOneContig);
                    } ).collect(Collectors.toList());
        }

        /**
         * Converts alignment records of the contig pointed to by {@code contigIdx} in a {@link FermiLiteAssembly} to custom {@link AlignmentInterval} format.
         * Filters out unmapped and ambiguous mappings (i.e. XA).
         */
        @VisibleForTesting
        private static List<AlignmentInterval> getAlignmentsForOneContig(final String contigName,
                                                                         final byte[] contigSequence,
                                                                         final List<BwaMemAlignment> contigAlignments,
                                                                         final List<String> refNames,
                                                                         final SAMFileHeader header) {

            return contigAlignments.stream()
                    .filter( bwaMemAlignment ->  bwaMemAlignment.getRefId() >= 0
                            && SAMFlag.SECONDARY_ALIGNMENT.isUnset(bwaMemAlignment.getSamFlag())) // mapped and not XA (i.e. not secondary)
                    .map(bwaMemAlignment -> BwaMemAlignmentUtils.applyAlignment(contigName, contigSequence, null,
                            null, bwaMemAlignment, refNames, header, false, false))
                    .map(AlignmentInterval::new)
                    .map(ar -> ContigAlignmentsModifier.splitGappedAlignment(ar, StructuralVariationDiscoveryArgumentCollection
                            .DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
                            .GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, contigSequence.length))
                    .flatMap(Utils::stream).collect(Collectors.toList());
        }

        @Override
        public JavaRDD<AlignedContig> getAlignedContigs() {

            // here we have two options, one is going through the route "BwaMemAlignment -> SAM -> GATKRead -> SAM -> AlignmentInterval"
            //                           which is the route if the discovery pipeline is run by "FindBreakpointEvidenceSpark -> write sam file -> load sam file -> DiscoverVariantsFromContigAlignmentsSAMSpark"
            //                         , the other is to go directly "BwaMemAlignment -> AlignmentInterval" by calling into {@code filterAndConvertToAlignedContigDirect()}, which is faster but not used here.
            //                         ; the two routes are tested to be generating the same output via {@code AlignedContigGeneratorUnitTest#testConvertAlignedAssemblyOrExcuseToAlignedContigsDirectAndConcordanceWithSAMRoute()}
            return filterAndConvertToAlignedContigViaSAM(alignedAssemblyOrExcuseList, header, ctx);
        }
    }

    public static Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls(final JavaSparkContext ctx,
                                                                              final SAMFileHeader header,
                                                                              final String cnvCallsFile) {
        final SVIntervalTree<VariantContext> cnvCalls;
        if (cnvCallsFile != null) {
            cnvCalls = CNVInputReader.loadCNVCalls(cnvCallsFile, header);
        } else {
            cnvCalls = null;
        }

        final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls;
        if (cnvCalls != null) {
            broadcastCNVCalls = ctx.broadcast(cnvCalls);
        } else {
            broadcastCNVCalls = null;
        }
        return broadcastCNVCalls;
    }

}
