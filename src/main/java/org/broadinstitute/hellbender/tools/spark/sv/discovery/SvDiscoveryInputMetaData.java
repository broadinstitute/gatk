package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;

import java.util.List;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;


public final class SvDiscoveryInputMetaData {

    public static final class ReferenceData {
        public final Broadcast<Set<String>> canonicalChromosomesBroadcast;
        public final Broadcast<ReferenceMultiSource> referenceBroadcast;
        public final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast;

        ReferenceData(final Broadcast<Set<String>> canonicalChromosomesBroadcast,
                      final Broadcast<ReferenceMultiSource> referenceBroadcast,
                      final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast) {
            this.canonicalChromosomesBroadcast = canonicalChromosomesBroadcast;
            this.referenceBroadcast = referenceBroadcast;
            this.referenceSequenceDictionaryBroadcast = referenceSequenceDictionaryBroadcast;
        }
    }

    public static final class SampleSpecificData {
        public final String sampleId;

        public final ReadMetadata readMetadata;
        public final Broadcast<SAMFileHeader> headerBroadcast;
        public final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast;
        public final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks;
        public final List<SVInterval> assembledIntervals;

        public SampleSpecificData(final String sampleId, final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast,
                                  final List<SVInterval> assembledIntervals,
                                  final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                  final ReadMetadata readMetadata,
                                  final Broadcast<SAMFileHeader> headerBroadcast) {
            this.sampleId = sampleId;
            this.cnvCallsBroadcast = cnvCallsBroadcast;
            this.assembledIntervals = assembledIntervals;
            this.evidenceTargetLinks = evidenceTargetLinks;
            this.readMetadata = readMetadata;
            this.headerBroadcast = headerBroadcast;
        }
    }

    public final ReferenceData referenceData;

    public final SampleSpecificData sampleSpecificData;

    public final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs;

    public final Logger toolLogger;

    public String outputPath;

    public SvDiscoveryInputMetaData(final JavaSparkContext ctx,
                                    final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs,
                                    final String nonCanonicalChromosomeNamesFile,
                                    final String outputPath,
                                    final ReadMetadata readMetadata,
                                    final List<SVInterval> assembledIntervals,
                                    final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                    final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast,
                                    final SAMFileHeader headerForReads,
                                    final ReferenceMultiSource reference,
                                    final Logger toolLogger) {

        final SAMSequenceDictionary sequenceDictionary = headerForReads.getSequenceDictionary();
        final Broadcast<Set<String>> canonicalChromosomesBroadcast =
                ctx.broadcast(SvDiscoveryUtils.getCanonicalChromosomes(nonCanonicalChromosomeNamesFile, sequenceDictionary));
        final String sampleId = SVUtils.getSampleId(headerForReads);

        this.referenceData = new ReferenceData(canonicalChromosomesBroadcast, ctx.broadcast(reference), ctx.broadcast(sequenceDictionary));
        this.sampleSpecificData = new SampleSpecificData(sampleId, cnvCallsBroadcast, assembledIntervals, evidenceTargetLinks, readMetadata, ctx.broadcast(headerForReads));
        this.discoverStageArgs = discoverStageArgs;
        this.outputPath = outputPath;
        this.toolLogger = toolLogger;
    }

    public void updateOutputPath(final String newOutputPath) {
        outputPath = newOutputPath;
    }

}
