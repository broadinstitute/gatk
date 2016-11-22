package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.collect.Iterables;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.Strandedness.SAME_STRAND;

@CommandLineProgramProperties(summary="Filter breakpoint alignments and call variants.",
        oneLineSummary="Filter breakpoint alignments and call variants",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class CallVariantsFromAlignedContigsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(CallVariantsFromAlignedContigsSpark.class);

    @Argument(doc = "URI of the output path", shortName = "outputPath",
            fullName = "outputPath", optional = false)
    private String outputPath;

    @Argument(doc = "output file name for called variants", shortName = "outputName",
            fullName = "outputName", optional = true)
    private String outputName = SVConstants.INS_DEL_OUTPUT_VCF;

    @Argument(doc = "Input file of contig alignments", shortName = "inputAlignments",
            fullName = "inputAlignments", optional = false)
    private String inputAlignments;

    @Argument(doc = "Input file of assembled contigs", shortName = "inputAssemblies",
            fullName = "inputAssemblies", optional = false)
    private String inputAssemblies;

    // This class requires a reference parameter in 2bit format (to broadcast) and a reference in FASTA format
    // (to get a good sequence dictionary).
    // todo: document this better
    @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
            fullName = "fastaReference", optional = false)
    private String fastaReference;


    private static final boolean DEBUG_STATS = true;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final PipelineOptions pipelineOptions = getAuthenticatedGCSOptions();
        final SAMSequenceDictionary referenceSequenceDictionary = new ReferenceMultiSource(pipelineOptions, fastaReference, ReferenceWindowFunctions.IDENTITY_FUNCTION).getReferenceSequenceDictionary(null);

        final JavaPairRDD<Iterable<AlignmentRegion>, byte[]> alignmentRegionsWithContigSequences = AssemblyAlignmentParser.prepAlignmentRegionsForCalling(ctx, inputAlignments, inputAssemblies).cache();

        if (DEBUG_STATS) debugStats(alignmentRegionsWithContigSequences, inputAlignments);

        final JavaRDD<VariantContext> variants = callVariantsFromAlignmentRegions(ctx.broadcast(getReference()), ctx.broadcast(referenceSequenceDictionary), alignmentRegionsWithContigSequences).cache();
        log.info("Called " + variants.count() + " variants");

        SVVCFWriter.writeVCF(pipelineOptions, outputPath, outputName, fastaReference, variants);
    }

    /**
     * This method processes an RDD containing alignment regions, scanning for chimeric alignments which match a set of filtering
     * criteria, and emitting a list of VariantContexts representing SVs for split alignments that pass.
     *
     * The input RDD is of the form:
     * Tuple2 of:
     *     {@code Iterable<AlignmentRegion>} AlignmentRegion objects representing all alignments for one contig
     *     A byte array with the sequence content of the contig
     * @param broadcastReference The broadcast handle to the reference (used to populate reference bases)
     * @param samSequenceDictionaryBroadcast sam sequence dictionary used for sorting alignment regions on difference chromosomes
     * @param alignmentRegionsWithContigSequences A data structure as described above, where a list of AlignmentRegions and the sequence of the contig are keyed by a tuple of Assembly ID and Contig ID
     */
    static JavaRDD<VariantContext> callVariantsFromAlignmentRegions(final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                    final Broadcast<SAMSequenceDictionary> samSequenceDictionaryBroadcast,
                                                                    final JavaPairRDD<Iterable<AlignmentRegion>, byte[]> alignmentRegionsWithContigSequences) {

        return alignmentRegionsWithContigSequences.filter(pair -> Iterables.size(pair._1())>1) // filter out any contigs that has less than two alignment records
                .flatMap(AssemblyAlignmentParser::getChimericAlignmentsFromAlignmentRegions)                                // 1. AR -> {CA}
                .mapToPair(ca -> new Tuple2<>(new BreakpointAllele(ca, samSequenceDictionaryBroadcast.getValue()), ca))     // 2. CA -> BA
                .filter(CallVariantsFromAlignedContigsSpark::filterOutBreakpointAlleleForNow)                               // development artifact (filter out variants not interested in for the moment)
                .groupByKey().map(tuple2 -> SVVariantConsensusCall.callVariantsFromConsensus(tuple2, broadcastReference));  // 3. {consensus BA} -> BA
    }

    private static void debugStats(JavaPairRDD<Iterable<AlignmentRegion>, byte[]> alignmentRegionsWithContigSequences, final String outPrefix) {
        log.info(alignmentRegionsWithContigSequences.count() + " contigs");
        final long noARs = alignmentRegionsWithContigSequences.filter(tuple2 -> Iterables.size(tuple2._1())==0).count();
        log.info(noARs + " contigs have no alignment regions");
        final long oneARs = alignmentRegionsWithContigSequences.filter(tuple2 -> Iterables.size(tuple2._1())==1).count();
        log.info(oneARs + " contigs have only one alignment regions");

        final JavaPairRDD<String, List<Tuple2<Integer, Integer>>> x = alignmentRegionsWithContigSequences.filter(tuple2 -> Iterables.size(tuple2._1())==2).mapToPair(tuple2 -> {
            final Iterator<AlignmentRegion> it = tuple2._1().iterator();
            final AlignmentRegion region1 = it.next(), region2 = it.next();
            return new Tuple2<>(region1.assemblyId+":"+region1.contigId, Arrays.asList(new Tuple2<>(region1.mapQual, region1.referenceInterval.size()), new Tuple2<>(region2.mapQual, region2.referenceInterval.size())));
        });
        x.coalesce(1).saveAsTextFile(outPrefix+"_2AR_MQs");
        log.info(x.count() + " contigs have two alignment regions");

        final JavaPairRDD<String, List<Tuple2<Integer, Integer>>> y = alignmentRegionsWithContigSequences.filter(tuple2 -> Iterables.size(tuple2._1())>2).mapToPair(tuple2 -> {
            final AlignmentRegion region1 = tuple2._1().iterator().next();
            return new Tuple2<>(region1.assemblyId+":"+region1.contigId, StreamSupport.stream(tuple2._1().spliterator(), false).map(ar -> new Tuple2<>(ar.mapQual, ar.referenceInterval.size())).collect(Collectors.toList()));
        });
        y.coalesce(1).saveAsTextFile(outPrefix+"_2+AR_MQs");
        log.info(y.count() + " contigs have more than two alignment regions");
    }

    // development artifact: filter out breakpoint alleles representing variants we don't care for now
    private static boolean filterOutBreakpointAlleleForNow(final Tuple2<BreakpointAllele, ChimericAlignment> breakpoint) {
        final BreakpointAllele breakpointAllele = breakpoint._1();
        final boolean isSameChromosome = breakpointAllele.leftJustifiedLeftBreakpoint.getContig().equals(breakpointAllele.leftJustifiedRightBreakpoint.getContig());
        final boolean noStrandSwitch = breakpointAllele.determineStrandedness() == SAME_STRAND;
        return isSameChromosome && noStrandSwitch;
    }
}
