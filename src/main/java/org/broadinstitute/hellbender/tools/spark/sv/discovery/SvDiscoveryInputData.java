package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;


public final class SvDiscoveryInputData {

    public final String sampleId;
    public String outputPath;
    public final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs;

    public final JavaRDD<GATKRead> assemblyRawAlignments;

    public final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast;
    public final List<SVInterval> assembledIntervals;
    public final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks;
    public final ReadMetadata metadata;

    public final Broadcast<SAMFileHeader> headerBroadcast;
    public final Broadcast<ReferenceMultiSource> referenceBroadcast;
    public final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast;

    public final Logger toolLogger;

    public SvDiscoveryInputData(final JavaSparkContext ctx,
                                final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs,
                                final String outputPath,
                                final ReadMetadata metadata,
                                final List<SVInterval> assembledIntervals,
                                final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast,
                                final JavaRDD<GATKRead> reads,
                                final SAMFileHeader headerForReads,
                                final ReferenceMultiSource reference,
                                final Logger toolLogger) {

        this(SVUtils.getSampleId(headerForReads), discoverStageArgs, outputPath, metadata, assembledIntervals,
                evidenceTargetLinks, reads, toolLogger,
                ctx.broadcast(reference), ctx.broadcast(headerForReads.getSequenceDictionary()), ctx.broadcast(headerForReads),
                cnvCallsBroadcast);
    }

    public SvDiscoveryInputData(final String sampleId,
                                final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs,
                                final String outputPath,
                                final ReadMetadata metadata,
                                final List<SVInterval> assembledIntervals,
                                final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                final JavaRDD<GATKRead> reads,
                                final Logger toolLogger,
                                final Broadcast<ReferenceMultiSource> referenceBroadcast,
                                final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast,
                                final Broadcast<SAMFileHeader> headerBroadcast,
                                final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast) {

        Utils.validate(! (evidenceTargetLinks != null && metadata == null),
                "Must supply read metadata when incorporating evidence target links");

        this.sampleId = sampleId;
        this.outputPath = outputPath;
        this.discoverStageArgs = discoverStageArgs;
        this.assemblyRawAlignments = reads;

        this.headerBroadcast = headerBroadcast;
        this.referenceBroadcast = referenceBroadcast;
        this.referenceSequenceDictionaryBroadcast = referenceSequenceDictionaryBroadcast;

        this.cnvCallsBroadcast = cnvCallsBroadcast;
        this.assembledIntervals = assembledIntervals;
        this.evidenceTargetLinks = evidenceTargetLinks;
        this.metadata = metadata;

        this.toolLogger = toolLogger;
    }

    public void updateOutputPath(final String newOutputPath) {
        outputPath = newOutputPath;
    }

}
