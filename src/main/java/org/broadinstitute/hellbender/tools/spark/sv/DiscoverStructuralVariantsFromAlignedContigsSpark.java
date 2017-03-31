package org.broadinstitute.hellbender.tools.spark.sv;

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

@CommandLineProgramProperties(summary="Filter breakpoint alignments and call variants.",
        oneLineSummary="Filter breakpoint alignments and call variants",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class DiscoverStructuralVariantsFromAlignedContigsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(DiscoverStructuralVariantsFromAlignedContigsSpark.class);

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
                = AssemblyAlignmentParser.prepAlignmentRegionsForCalling(ctx, inputAlignments, inputAssemblies).cache();

        final Broadcast<ReferenceMultiSource> broadcastReference = ctx.broadcast(getReference());

        final JavaRDD<VariantContext> variants
                = SVVariantConsensusDiscovery.discoverNovelAdjacencyFromChimericAlignments(alignmentRegionsWithContigSequences, log)
                .map(tuple2 -> SVVariantConsensusDiscovery.discoverVariantsFromConsensus(tuple2, broadcastReference));

        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputPath, outputName, fastaReference, variants, log);
    }

}
