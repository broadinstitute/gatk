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
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

/**
 * This tool takes a SAM file containing the alignments of assembled contigs or long reads to the reference
 * and searches it for split alignments indicating the presence of structural variations. To do so the tool parses
 * primary and supplementary alignments; secondary alignments are ignored. To be considered valid evidence of an SV,
 * two alignments from the same contig must have mapping quality 60, and both alignments must have length greater than
 * or equal to minAlignmentLength.
 */
@CommandLineProgramProperties(summary="Parse a SAM file containing contigs or long reads aligned to the reference, and call SVs",
        oneLineSummary="Parse a SAM file containing contigs or long reads aligned to the reference, and call SVs",
        programGroup = SparkProgramGroup.class)
public final class DiscoverStructuralVariantsFromAlignedContigsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(DiscoverStructuralVariantsFromAlignedContigsSAMSpark.class);

    @Argument(doc = "URL of the output path", shortName = "outputPath",
            fullName = "outputPath", optional = false)
    private String outputPath;

    @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
            fullName = "fastaReference", optional = false)
    private String fastaReference;

    @Argument(doc = "Minimum flanking alignment length", shortName = "minAlignLength",
            fullName = "minAlignLength", optional = true)
    private Integer minAlignLength = SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH;

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

        final JavaPairRDD<Iterable<AlignmentRegion>, byte[]> alignmentRegionsIterable
                = getReads().groupBy(GATKRead::getName).map(Tuple2::_2).mapToPair(AssemblyAlignmentParser::convertToAlignmentRegions);

        final Broadcast<ReferenceMultiSource> broadcastReference = ctx.broadcast(getReference());

        final JavaRDD<VariantContext> variants
                = SVVariantConsensusDiscovery.discoverNovelAdjacencyFromChimericAlignments(alignmentRegionsIterable, log)
                .map(tuple2 -> SVVariantConsensusDiscovery.discoverVariantsFromConsensus(tuple2, broadcastReference));

        SVVCFWriter.writeVCF(getAuthenticatedGCSOptions(), outputPath, SVConstants.DiscoveryStepConstants.CURRENTLY_CAPABLE_VARIANTS_VCF, fastaReference, variants, log);
    }

}
