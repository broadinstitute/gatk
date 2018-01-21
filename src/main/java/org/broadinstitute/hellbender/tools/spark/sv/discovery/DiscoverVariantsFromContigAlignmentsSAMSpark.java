package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
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
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyReferenceLocations;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleNovelAdjacencyInterpreter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import scala.Tuple2;

import java.util.List;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH;

/**
 * (Internal) Examines aligned contigs from local assemblies and calls structural variants
 *
 * <p>This tool is used in development and should not be of interest to most researchers.  It packages structural
 * variant calling as a separate tool, independent of the generation of local assemblies.
 * Most researchers will run StructuralVariationDiscoveryPipelineSpark, which both generates local assemblies
 * of interesting genomic regions, and then calls structural variants from these assemblies.</p>
 * <p>This tool takes a SAM/BAM/CRAM containing the alignments of assembled contigs from local assemblies
 * and searches it for split alignments indicating the presence of structural variations. To do so the tool parses
 * primary and supplementary alignments; secondary alignments are ignored. To be considered valid evidence of an SV,
 * two alignments from the same contig must have mapping quality 60, and both alignments must have length greater than
 * or equal to min-alignment-length. Imprecise variants with approximate locations are also called.</p>
 * <p>The input file is typically the output file produced by FindBreakpointEvidenceSpark.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of assembled contigs or long reads aligned to reference.</li>
 *     <li>The reference to which the contigs have been aligned.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A vcf file describing the discovered structural variants.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk DiscoverVariantsFromContigAlignmentsSAMSpark \
 *     -I assemblies.sam \
 *     -R reference.2bit \
 *     -O structural_variants.vcf
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 *
 * <h3>Notes</h3>
 * <p>The reference is broadcast by Spark, and must therefore be a .2bit file due to current restrictions.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Examines aligned contigs from local assemblies and calls structural variants",
        summary =
        "This tool is used in development and should not be of interest to most researchers.  It packages structural" +
        " variant calling as a separate tool, independent of the generation of local assemblies." +
        " Most researchers will run StructuralVariationDiscoveryPipelineSpark, which both generates local assemblies" +
        " of interesting genomic regions, and then calls structural variants from these assemblies." +
        " This tool takes a SAM/BAM/CRAM containing the alignments of assembled contigs from local assemblies" +
        " and searches it for split alignments indicating the presence of structural variations. To do so the tool parses" +
        " primary and supplementary alignments; secondary alignments are ignored. To be considered valid evidence of an SV," +
        " two alignments from the same contig must have mapping quality 60, and both alignments must have length greater than" +
        " or equal to min-alignment-length. Imprecise variants with approximate locations are also called.\n" +
        " The input file is typically the output file produced by FindBreakpointEvidenceSpark.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class DiscoverVariantsFromContigAlignmentsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(DiscoverVariantsFromContigAlignmentsSAMSpark.class);

    @ArgumentCollection
    private final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs =
            new DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

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

        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast =
                StructuralVariationDiscoveryPipelineSpark.broadcastCNVCalls(ctx, getHeaderForReads(),
                        discoverStageArgs.cnvCallsFile);

        final SvDiscoveryInputData svDiscoveryInputData =
                new SvDiscoveryInputData(ctx, discoverStageArgs, vcfOutputFileName,
                        null, null, null,
                        cnvCallsBroadcast,
                        getReads(), getHeaderForReads(), getReference(), localLogger);

        final JavaRDD<AlignedContig> parsedContigAlignments =
                new SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SAMFormattedContigAlignmentParser(svDiscoveryInputData.assemblyRawAlignments, svDiscoveryInputData.headerBroadcast.getValue(), true)
                        .getAlignedContigs();

        // assembly-based breakpoints
        List<VariantContext> annotatedVariants = discoverVariantsFromChimeras(svDiscoveryInputData, parsedContigAlignments);

        final SAMSequenceDictionary refSeqDictionary = svDiscoveryInputData.referenceSequenceDictionaryBroadcast.getValue();
        SVVCFWriter.writeVCF(annotatedVariants, vcfOutputFileName, refSeqDictionary, localLogger);
    }

    @Deprecated
    public static List<VariantContext> discoverVariantsFromChimeras(final SvDiscoveryInputData svDiscoveryInputData,
                                                                    final JavaRDD<AlignedContig> alignedContigs) {
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;

        final JavaPairRDD<byte[], List<ChimericAlignment>> contigSeqAndChimeras =
                alignedContigs.filter(alignedContig -> alignedContig.alignmentIntervals.size() > 1)
                        .mapToPair(alignedContig ->
                                new Tuple2<>(alignedContig.contigSequence,
                                        ChimericAlignment.parseOneContig(alignedContig, referenceSequenceDictionaryBroadcast.getValue(),
                                                true, DEFAULT_MIN_ALIGNMENT_LENGTH,
                                                CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, true)));

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;
        final List<SVInterval> assembledIntervals = svDiscoveryInputData.assembledIntervals;
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputData.cnvCallsBroadcast;
        final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputData.discoverStageArgs;
        final String sampleId = svDiscoveryInputData.sampleId;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        final JavaPairRDD<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> narlsAndSources =
                contigSeqAndChimeras
                        .flatMapToPair(tigSeqAndChimeras -> {
                            final byte[] contigSeq = tigSeqAndChimeras._1;
                            final List<ChimericAlignment> chimericAlignments = tigSeqAndChimeras._2;
                            final Stream<Tuple2<NovelAdjacencyReferenceLocations, ChimericAlignment>> novelAdjacencyAndSourceChimera =
                                    chimericAlignments.stream()
                                            .map(ca -> new Tuple2<>(
                                                    new NovelAdjacencyReferenceLocations(ca, contigSeq,
                                                            referenceSequenceDictionaryBroadcast.getValue()), ca));
                            return novelAdjacencyAndSourceChimera.iterator();
                        })
                        .groupByKey()   // group the same novel adjacency produced by different contigs together
                        .cache();


        try {
            SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals, narlsAndSources.map(Tuple2::_1).collect(),
                    referenceSequenceDictionaryBroadcast.getValue(), discoverStageArgs, toolLogger);

            List<VariantContext> annotatedVariants =
                    narlsAndSources
                            .mapToPair(noveltyAndEvidence -> new Tuple2<>(noveltyAndEvidence._1,
                                    new Tuple2<>(SimpleNovelAdjacencyInterpreter.inferSimpleTypeFromNovelAdjacency(noveltyAndEvidence._1), noveltyAndEvidence._2)))       // type inference based on novel adjacency and evidence alignments
                            .map(noveltyTypeAndEvidence ->
                            {
                                final NovelAdjacencyReferenceLocations novelAdjacency = noveltyTypeAndEvidence._1;
                                final SimpleSVType inferredSimpleType = noveltyTypeAndEvidence._2._1;
                                final Iterable<ChimericAlignment> evidence = noveltyTypeAndEvidence._2._2;
                                return AnnotatedVariantProducer
                                        .produceAnnotatedVcFromInferredTypeAndRefLocations(
                                                novelAdjacency.leftJustifiedLeftRefLoc,
                                                novelAdjacency.leftJustifiedRightRefLoc.getStart(),
                                                novelAdjacency.complication,
                                                inferredSimpleType,
                                                null, // TODO: 1/21/18 implement this for InsDel
                                                evidence,
                                                referenceBroadcast,
                                                referenceSequenceDictionaryBroadcast,
                                                cnvCallsBroadcast,
                                                sampleId);
                            })
                            .collect();
            return annotatedVariants;
        } finally {
            narlsAndSources.unpersist();
        }
    }

}
