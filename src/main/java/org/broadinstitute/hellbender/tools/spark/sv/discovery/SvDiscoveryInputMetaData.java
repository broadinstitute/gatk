package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
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
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.List;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection;


public final class SvDiscoveryInputMetaData {

    public static final class ReferenceData {
        private final Broadcast<Set<String>> canonicalChromosomesBroadcast;
        private final Broadcast<ReferenceMultiSource> referenceBroadcast;
        private final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast;

        // marked private so can be constructed only by enclosing class
        private ReferenceData(@Nonnull final Broadcast<Set<String>> canonicalChromosomesBroadcast,
                              @Nonnull final Broadcast<ReferenceMultiSource> referenceBroadcast,
                              @Nonnull final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast) {
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

        private final ReadMetadata readMetadata;
        private final Broadcast<SAMFileHeader> headerBroadcast;
        private final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast;
        private final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks;
        private final List<SVInterval> assembledIntervals;

        // marked private so can be constructed only by enclosing class

        /**
         * Holds sample specific data.
         * @param sampleId              sample ID string
         * @param headerBroadcast       SAM file header corresponding to the sample
         * @param cnvCallsBroadcast     broadcast of CNV calls from CNV caller(s), can be null
         * @param assembledIntervals    broadcast of assembled intervals, can be null
         * @param evidenceTargetLinks   broadcast of evidence target links, can be null
         * @param readMetadata          broadcast of short read meta data, can be null
         */
        private SampleSpecificData(@Nonnull final String sampleId, @Nonnull final Broadcast<SAMFileHeader> headerBroadcast,
                                   @Nullable final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast,
                                   @Nullable final List<SVInterval> assembledIntervals,
                                   @Nullable final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                   @Nullable final ReadMetadata readMetadata) {
            this.sampleId = sampleId;
            this.cnvCallsBroadcast = cnvCallsBroadcast;
            this.assembledIntervals = assembledIntervals;
            this.evidenceTargetLinks = evidenceTargetLinks;
            this.readMetadata = readMetadata;
            this.headerBroadcast = headerBroadcast;
        }

        public String getSampleId() {
            return sampleId;
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

    private final Set<VCFHeaderLine> defaultToolVCFHeaderLines;

    private final Logger toolLogger;

    private String outputPath;

    /**
     * Holds various information used in the discovery stage.
     * @param ctx                               for broadcasting
     * @param discoverStageArgs                 argument collection used for discovery
     * @param reference                         reference bases
     * @param outputPath                        path to output files for discovery stage
     * @param headerForReads                    SAM file header for the sample's bam
     * @param defaultToolVCFHeaderLines         a set of default VCF headers to be included in the output VCF
     * @param toolLogger                        logger
     * @param nonCanonicalChromosomeNamesFile   path to file holding non-canonical chromosome names, all chromosomes considered to be canonical when {@code null}
     * @param readMetadata                      short read meta data, can be null
     * @param assembledIntervals                assembled intervals, can be null
     * @param evidenceTargetLinks               evidence target links, can be null
     * @param cnvCallsBroadcast                 CNV calls by CNV caller(s), can be null
     */
    public SvDiscoveryInputMetaData(@Nonnull final JavaSparkContext ctx,
                                    @Nonnull final DiscoverVariantsFromContigAlignmentsSparkArgumentCollection discoverStageArgs,
                                    @Nonnull final ReferenceMultiSource reference,
                                    @Nonnull final String outputPath,
                                    @Nonnull final SAMFileHeader headerForReads,
                                    @Nonnull final Set<VCFHeaderLine> defaultToolVCFHeaderLines,
                                    @Nonnull final Logger toolLogger,
                                    @Nullable final String nonCanonicalChromosomeNamesFile,
                                    @Nullable final ReadMetadata readMetadata,
                                    @Nullable final List<SVInterval> assembledIntervals,
                                    @Nullable final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                    @Nullable final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast) {

        Utils.nonNull(ctx);
        Utils.nonNull(discoverStageArgs);
        Utils.nonNull(headerForReads);
        Utils.nonNull(reference);
        Utils.nonNull(defaultToolVCFHeaderLines);
        Utils.nonNull(outputPath);
        Utils.nonNull(toolLogger);

        final SAMSequenceDictionary sequenceDictionary = headerForReads.getSequenceDictionary();
        final Broadcast<Set<String>> canonicalChromosomesBroadcast =
                ctx.broadcast(SVUtils.getCanonicalChromosomes(nonCanonicalChromosomeNamesFile, sequenceDictionary));
        final String sampleId = SVUtils.getSampleId(headerForReads);

        this.referenceData = new ReferenceData(canonicalChromosomesBroadcast, ctx.broadcast(reference), ctx.broadcast(sequenceDictionary));
        this.sampleSpecificData = new SampleSpecificData(sampleId, ctx.broadcast(headerForReads), cnvCallsBroadcast, assembledIntervals, evidenceTargetLinks, readMetadata);
        this.discoverStageArgs = discoverStageArgs;
        this.outputPath = outputPath;
        this.defaultToolVCFHeaderLines = defaultToolVCFHeaderLines;
        this.toolLogger = toolLogger;
    }


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

    public Set<VCFHeaderLine> getDefaultToolVCFHeaderLines() {
        return defaultToolVCFHeaderLines;
    }

    public void updateOutputPath(final String newOutputPath) {
        outputPath = newOutputPath;
    }

}
