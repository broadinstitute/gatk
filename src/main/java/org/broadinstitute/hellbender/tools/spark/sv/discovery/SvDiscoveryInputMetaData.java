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

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection;


public final class SvDiscoveryInputMetaData {

    public ReferenceData getReferenceData() {
        return referenceData;
    }

    public SampleSpecificData getSampleSpecificData() {
        return sampleSpecificData;
    }

    public DiscoverVariantsFromContigAlignmentsSparkArgumentCollection getDiscoverStageArgs() {
        return discoverStageArgs;
    }

    public Logger getToolLogger() {
        return toolLogger;
    }

    public String getOutputPath() {
        return outputPath;
    }

    public static final class ReferenceData {
        private final Broadcast<Set<String>> canonicalChromosomesBroadcast;
        private final Broadcast<ReferenceMultiSource> referenceBroadcast;
        private final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast;

        ReferenceData(final Broadcast<Set<String>> canonicalChromosomesBroadcast,
                      final Broadcast<ReferenceMultiSource> referenceBroadcast,
                      final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast) {
            this.canonicalChromosomesBroadcast = canonicalChromosomesBroadcast;
            this.referenceBroadcast = referenceBroadcast;
            this.referenceSequenceDictionaryBroadcast = referenceSequenceDictionaryBroadcast;
        }

        public Broadcast<Set<String>> getCanonicalChromosomesBroadcast() {
            return canonicalChromosomesBroadcast;
        }

        public Broadcast<ReferenceMultiSource> getReferenceBroadcast() {
            return referenceBroadcast;
        }

        public Broadcast<SAMSequenceDictionary> getReferenceSequenceDictionaryBroadcast() {
            return referenceSequenceDictionaryBroadcast;
        }
    }

    public static final class SampleSpecificData {
        private final String sampleId;

        private final String pathToAssemblyFastqDir;

        private final ReadMetadata readMetadata;
        private final Broadcast<SAMFileHeader> headerBroadcast;
        private final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast;
        private final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks;
        private final List<SVInterval> assembledIntervals;

        public SampleSpecificData(final String sampleId, final String pathToAssemblyFastqDir,
                                  final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast,
                                  final List<SVInterval> assembledIntervals,
                                  final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                  final ReadMetadata readMetadata,
                                  final Broadcast<SAMFileHeader> headerBroadcast) {
            this.sampleId = sampleId;
            this.pathToAssemblyFastqDir = pathToAssemblyFastqDir;
            this.cnvCallsBroadcast = cnvCallsBroadcast;
            this.assembledIntervals = assembledIntervals;
            this.evidenceTargetLinks = evidenceTargetLinks;
            this.readMetadata = readMetadata;
            this.headerBroadcast = headerBroadcast;
        }

        public String getSampleId() {
            return sampleId;
        }

        public String getPathToAssemblyFastqDir() {
            return pathToAssemblyFastqDir;
        }

        public ReadMetadata getReadMetadata() {
            return readMetadata;
        }

        public Broadcast<SAMFileHeader> getHeaderBroadcast() {
            return headerBroadcast;
        }

        public Broadcast<SVIntervalTree<VariantContext>> getCnvCallsBroadcast() {
            return cnvCallsBroadcast;
        }

        public PairedStrandedIntervalTree<EvidenceTargetLink> getEvidenceTargetLinks() {
            return evidenceTargetLinks;
        }

        public List<SVInterval> getAssembledIntervals() {
            return assembledIntervals;
        }
    }

    private final ReferenceData referenceData;

    private final SampleSpecificData sampleSpecificData;

    private final DiscoverVariantsFromContigAlignmentsSparkArgumentCollection discoverStageArgs;

    private final Logger toolLogger;

    private String outputPath;

    public SvDiscoveryInputMetaData(final JavaSparkContext ctx,
                                    final DiscoverVariantsFromContigAlignmentsSparkArgumentCollection discoverStageArgs,
                                    final String nonCanonicalChromosomeNamesFile,
                                    final String outputPath,
                                    final ReadMetadata readMetadata,
                                    final List<SVInterval> assembledIntervals,
                                    final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                    final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast,
                                    final String pathToAssemblyFastqDir,
                                    final SAMFileHeader headerForReads,
                                    final ReferenceMultiSource reference,
                                    final Logger toolLogger) {

        final SAMSequenceDictionary sequenceDictionary = headerForReads.getSequenceDictionary();
        final Broadcast<Set<String>> canonicalChromosomesBroadcast =
                ctx.broadcast(SvDiscoveryUtils.getCanonicalChromosomes(nonCanonicalChromosomeNamesFile, sequenceDictionary));
        final String sampleId = SVUtils.getSampleId(headerForReads);

        this.referenceData = new ReferenceData(canonicalChromosomesBroadcast, ctx.broadcast(reference), ctx.broadcast(sequenceDictionary));
        this.sampleSpecificData = new SampleSpecificData(sampleId, pathToAssemblyFastqDir, cnvCallsBroadcast, assembledIntervals, evidenceTargetLinks, readMetadata, ctx.broadcast(headerForReads));
        this.discoverStageArgs = discoverStageArgs;
        this.outputPath = outputPath;
        this.toolLogger = toolLogger;
    }

    public void updateOutputPath(final String newOutputPath) {
        outputPath = newOutputPath;
    }

}
