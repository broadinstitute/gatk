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
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.AssemblyContigAlignmentSignatureClassifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.FindBreakpointEvidenceSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Serializable;

import java.util.EnumMap;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;

/**
 * Runs the structural variation discovery workflow on a single sample
 *
 * <p>This tool packages the algorithms described in {@link FindBreakpointEvidenceSpark} and
 * {@link DiscoverVariantsFromContigAlignmentsSAMSpark} as an integrated workflow.  Please consult the
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

    @Argument(doc = "sam file for aligned contigs", fullName = "contig-sam-file")
    private String outputAssemblyAlignments;

    @Argument(doc = "filename for output vcf", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutputFileName;

    @ArgumentCollection
    private final FindBreakpointEvidenceSparkArgumentCollection evidenceAndAssemblyArgs
            = new FindBreakpointEvidenceSparkArgumentCollection();

    @ArgumentCollection
    private final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs
            = new DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Advanced
    @Argument(doc = "directory to output results of our prototyping breakpoint and type inference tool in addition to the master VCF;" +
            " the directory contains multiple VCF's for different types and record-generating SAM files of assembly contigs,",
            fullName = "exp-variants-out-dir", optional = true)
    private String expVariantsOutDir;

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
        final Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(headerForReads);
        final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast = ctx.broadcast(getReference());
        final SAMSequenceDictionary refSequenceDictionary = headerForReads.getSequenceDictionary();
        final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary = ctx.broadcast(refSequenceDictionary);
        final String sampleId = SVUtils.getSampleId(headerForReads);

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

        final List<EvidenceTargetLink> evidenceTargetLinks = assembledEvidenceResults.getEvidenceTargetLinks();
        final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceLinkTree = makeEvidenceLinkTree(evidenceTargetLinks);

        final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls =
                DiscoverVariantsFromContigAlignmentsSAMSpark.broadcastCNVCalls(ctx, headerForReads, sampleId, discoverStageArgs);

        // discover variants and write to vcf
        DiscoverVariantsFromContigAlignmentsSAMSpark
                .discoverVariantsAndWriteVCF(
                        parsedAlignments,
                        assembledEvidenceResults.getAssembledIntervals(),
                        discoverStageArgs,
                        referenceMultiSourceBroadcast,
                        broadcastSequenceDictionary,
                        vcfOutputFileName,
                        localLogger,
                        evidenceLinkTree,
                        assembledEvidenceResults.getReadMetadata(),
                        broadcastCNVCalls,
                        sampleId
                );

        if ( expVariantsOutDir != null ) {

            experimentalInterpretation(ctx,
                    headerForReads,
                    headerBroadcast,
                    referenceMultiSourceBroadcast,
                    refSequenceDictionary,
                    broadcastSequenceDictionary,
                    sampleId,
                    assembledEvidenceResults);
        }
    }

    // hook up prototyping breakpoint and type inference tool
    private void experimentalInterpretation(final JavaSparkContext ctx,
                                            final SAMFileHeader headerForReads,
                                            final Broadcast<SAMFileHeader> headerBroadcast,
                                            final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast,
                                            final SAMSequenceDictionary refSequenceDictionary,
                                            final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                            final String sampleId,
                                            final FindBreakpointEvidenceSpark.AssembledEvidenceResults assembledEvidenceResults) {

        if ( expVariantsOutDir == null )
            return;

        final SAMReadGroupRecord contigAlignmentsReadGroup = new SAMReadGroupRecord(SVUtils.GATKSV_CONTIG_ALIGNMENTS_READ_GROUP_ID);
        final List<String> refNames = SequenceDictionaryUtils.getContigNamesList(refSequenceDictionary);

        List<GATKRead> readsList =
                assembledEvidenceResults
                        .getAlignedAssemblyOrExcuseList().stream()
                        .filter(AlignedAssemblyOrExcuse::isNotFailure)
                        .flatMap(aa -> aa.toSAMStreamForAlignmentsOfThisAssembly(headerForReads, refNames, contigAlignmentsReadGroup))
                        .map(SAMRecordToGATKReadAdapter::new)
                        .collect(Collectors.toList());
        JavaRDD<GATKRead> reads = ctx.parallelize(readsList);

        EnumMap<AssemblyContigAlignmentSignatureClassifier.RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>>
                contigsByPossibleRawTypes =
                SvDiscoverFromLocalAssemblyContigAlignmentsSpark.preprocess(reads, headerBroadcast,
                broadcastSequenceDictionary, evidenceAndAssemblyArgs.crossContigsToIgnoreFile,
                        expVariantsOutDir, true, localLogger);

        SvDiscoverFromLocalAssemblyContigAlignmentsSpark.dispatchJobs(sampleId, expVariantsOutDir,
                contigsByPossibleRawTypes, referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger);

        // TODO: 11/30/17 add EXTERNAL_CNV_CALLS annotation to the variants called here
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

        @Override
        public JavaRDD<AlignedContig> getAlignedContigs() {

            // here we have two options, one is going through the route "BwaMemAlignment -> SAM -> GATKRead -> SAM -> AlignmentInterval"
            //                           which is the route if the discovery pipeline is run by "FindBreakpointEvidenceSpark -> write sam file -> load sam file -> DiscoverVariantsFromContigAlignmentsSAMSpark"
            //                         , the other is to go directly "BwaMemAlignment -> AlignmentInterval" by calling into {@code filterAndConvertToAlignedContigDirect()}, which is faster but not used here.
            //                         ; the two routes are tested to be generating the same output via {@code AlignedContigGeneratorUnitTest#testConvertAlignedAssemblyOrExcuseToAlignedContigsDirectAndConcordanceWithSAMRoute()}
            return filterAndConvertToAlignedContigViaSAM(alignedAssemblyOrExcuseList, header, ctx);
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
                            DiscoverVariantsFromContigAlignmentsSAMSpark.
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
                    .filter(contig -> !contig.alignmentIntervals.isEmpty())     // filter out unmapped and contigs without primary alignments
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
                        final List<AlignmentInterval> arOfAContig
                                = getAlignmentsForOneContig(contigName, contigSequence, allAlignments.get(contigIdx), refNames, header);
                        return new AlignedContig(contigName, contigSequence, arOfAContig, false);
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
                    .map(ar -> GappedAlignmentSplitter.split(ar, StructuralVariationDiscoveryArgumentCollection
                            .DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
                            .GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, contigSequence.length))
                    .flatMap(Utils::stream).collect(Collectors.toList());
        }
    }

}
