package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.collect.Iterables;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.spark.sv.ContigsCollection.loadContigsCollectionKeyedByAssemblyId;

@CommandLineProgramProperties(summary="Filter breakpoint alignments and call variants.",
        oneLineSummary="Filter breakpoint alignments and call variants",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class DiscoverVariantsFromAlignedSGAContigsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(DiscoverVariantsFromAlignedSGAContigsSpark.class);

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection  discoverVariantsFromContigsAlignmentsSparkArgumentCollection
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "URI of the output path", shortName = "outputPath",
            fullName = "outputPath", optional = false)
    private String outputPath;

    @Argument(doc = "output file name for called variants", shortName = "outputName",
            fullName = "outputName", optional = true)
    private String outputName = SVConstants.DiscoveryStepConstants.CURRENTLY_CAPABLE_VARIANTS_VCF;

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


    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaPairRDD<Iterable<AlignmentRegion>, byte[]> alignmentRegionsWithContigSequences
                = prepAlignmentRegionsForCalling(ctx, inputAlignments, inputAssemblies, discoverVariantsFromContigsAlignmentsSparkArgumentCollection.logContigAlignmentSimpleStats).cache();

        final Broadcast<ReferenceMultiSource> broadcastReference = ctx.broadcast(getReference());

        final JavaRDD<VariantContext> variants
                = SVVariantConsensusDiscovery.discoverNovelAdjacencyFromChimericAlignments(alignmentRegionsWithContigSequences, log)
                .map(tuple2 -> SVVariantConsensusDiscovery.discoverVariantsFromConsensus(tuple2, broadcastReference));

        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputName, fastaReference, variants, log);
    }

    /**
     * Loads the alignment regions and sequence of all locally-assembled contigs by SGA from the text file they are in;
     * one record for each contig.
     * @param ctx                       spark context for IO operations
     * @param pathToInputAlignments     path string to alignments of the contigs; format assumed to be consistent/parsable by {@link AlignmentRegion#toString()}
     * @param logContigAlignmentSimpleStats
     * @return                          an PairRDD for all assembled contigs with their alignment regions and sequence
     */
    private static JavaPairRDD<Iterable<AlignmentRegion>, byte[]> prepAlignmentRegionsForCalling(final JavaSparkContext ctx,
                                                                                                 final String pathToInputAlignments,
                                                                                                 final String pathToInputAssemblies,
                                                                                                 final boolean logContigAlignmentSimpleStats) {

        final JavaPairRDD<Tuple2<String, String>, Iterable<AlignmentRegion>> alignmentRegionsKeyedByAssemblyAndContigId
                = parseAlignments(ctx, pathToInputAlignments);
        if (logContigAlignmentSimpleStats) {
            debugStats(alignmentRegionsKeyedByAssemblyAndContigId, pathToInputAlignments);
        }

        final JavaPairRDD<Tuple2<String, String>, byte[]> contigSequences
                = loadContigSequence(ctx, pathToInputAssemblies);

        return alignmentRegionsKeyedByAssemblyAndContigId
                .join(contigSequences)
                .mapToPair(Tuple2::_2);
    }

    // TODO: 11/23/16 test
    private static JavaPairRDD<Tuple2<String, String>, Iterable<AlignmentRegion>> parseAlignments(final JavaSparkContext ctx,
                                                                                                  final String pathToInputAlignments) {

        return ctx.textFile(pathToInputAlignments).map(textLine -> AlignmentRegion.fromString(textLine.split(AlignmentRegion.STRING_REP_SEPARATOR,-1)))
                .flatMap(oneRegion -> AssemblyAlignmentParser.breakGappedAlignment(oneRegion, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY).iterator())
                .mapToPair(alignmentRegion -> new Tuple2<>(new Tuple2<>(alignmentRegion.assemblyId, alignmentRegion.contigId), alignmentRegion))
                .groupByKey();
    }

    // TODO: 12/15/16 test
    private static JavaPairRDD<Tuple2<String, String>, byte[]> loadContigSequence(final JavaSparkContext ctx,
                                                                                  final String pathToInputAssemblies) {
        return loadContigsCollectionKeyedByAssemblyId(ctx, pathToInputAssemblies)
                .flatMapToPair(assemblyIdAndContigsCollection -> {
                    final String assemblyId = assemblyIdAndContigsCollection._1;
                    final ContigsCollection contigsCollection = assemblyIdAndContigsCollection._2;
                    return contigsCollection.getContents().stream()
                            .map(pair -> new Tuple2<>(new Tuple2<>(assemblyId, pair._1.toString()), pair._2.toString().getBytes()))
                            .collect(Collectors.toList()).iterator();
                });
    }


    private static void debugStats(final JavaPairRDD<Tuple2<String, String>, Iterable<AlignmentRegion>> alignmentRegionsWithContigSequences, final String outPrefix) {
        log.info(alignmentRegionsWithContigSequences.count() + " contigs");
        final long noARs = alignmentRegionsWithContigSequences.filter(tuple2 -> Iterables.size(tuple2._2())==0).count();
        log.info(noARs + " contigs have no alignments");
        final long oneARs = alignmentRegionsWithContigSequences.filter(tuple2 -> Iterables.size(tuple2._2())==1).count();
        log.info(oneARs + " contigs have only one alignments");

        final JavaPairRDD<String, List<Tuple2<Integer, Integer>>> x = alignmentRegionsWithContigSequences.filter(tuple2 -> Iterables.size(tuple2._2())==2).mapToPair(tuple2 -> {
            final Iterator<AlignmentRegion> it = tuple2._2().iterator();
            final AlignmentRegion region1 = it.next(), region2 = it.next();
            return new Tuple2<>(region1.assemblyId+":"+region1.contigId, Arrays.asList(new Tuple2<>(region1.mapQual, region1.referenceInterval.size()), new Tuple2<>(region2.mapQual, region2.referenceInterval.size())));
        });
        x.coalesce(1).saveAsTextFile(outPrefix+"_withTwoAlignments");
        log.info(x.count() + " contigs have two alignments");

        final JavaPairRDD<String, List<Tuple2<Integer, Integer>>> y = alignmentRegionsWithContigSequences.filter(tuple2 -> Iterables.size(tuple2._2())>2).mapToPair(tuple2 -> {
            final AlignmentRegion region1 = tuple2._2().iterator().next();
            return new Tuple2<>(region1.assemblyId+":"+region1.contigId, StreamSupport.stream(tuple2._2().spliterator(), false).map(ar -> new Tuple2<>(ar.mapQual, ar.referenceInterval.size())).collect(Collectors.toList()));
        });
        y.coalesce(1).saveAsTextFile(outPrefix+"_withMoreThanTwoAlignments");
        log.info(y.count() + " contigs have more than two alignments");
    }

}
