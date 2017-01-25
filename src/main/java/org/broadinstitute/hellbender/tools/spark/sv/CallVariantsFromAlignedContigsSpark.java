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
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import scala.Tuple2;

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
    private String outputName = SVConstants.CallingStepConstants.CURRENTLY_CAPABLE_VARIANTS_VCF;

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
                = AssemblyAlignmentParser.prepAlignmentRegionsForCalling(ctx, inputAlignments, inputAssemblies).cache();

        final JavaRDD<VariantContext> variants
                = callVariantsFromAlignmentRegions(ctx.broadcast(getReference()), alignmentRegionsWithContigSequences, log).cache();
        log.info("Called " + variants.count() + " variants");

        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputPath, outputName, fastaReference, variants);
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
     * @param alignmentRegionsWithContigSequence A data structure as described above, where a list of AlignmentRegions and the sequence of the contig are keyed by a tuple of Assembly ID and Contig ID
     * @param logger
     */
    static JavaRDD<VariantContext> callVariantsFromAlignmentRegions(final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                    final JavaPairRDD<Iterable<AlignmentRegion>, byte[]> alignmentRegionsWithContigSequence,
                                                                    final Logger logger) {

        return alignmentRegionsWithContigSequence.filter(pair -> Iterables.size(pair._1())>1) // filter out any contigs that has less than two alignment records
                .flatMap( input -> ChimericAlignment.fromSplitAlignments(input).iterator())                                                            // 1. AR -> {CA}
                .mapToPair(ca -> new Tuple2<>(new NovelAdjacencyReferenceLocations(ca), ca))                                // 2. CA -> BP
                .groupByKey()                                                                                               // 3. {consensus BP}
                .map(tuple2 -> SVVariantConsensusCall.callVariantsFromConsensus(tuple2, broadcastReference));               // BP annotated with list of CA's
    }
}
