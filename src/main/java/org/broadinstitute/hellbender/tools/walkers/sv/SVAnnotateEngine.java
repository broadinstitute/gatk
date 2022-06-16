package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.compress.utils.Sets;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfTranscriptFeature;

import java.util.*;

public class SVAnnotateEngine {
    private final int maxBreakendLen;
    private final GTFIntervalTreesContainer gtfIntervalTrees;
    private final SVIntervalTree<String> nonCodingIntervalTree;
    private final SAMSequenceDictionary sequenceDictionary;

    private final Set<String> MSV_EXON_OVERLAP_CLASSIFICATIONS = Sets.newHashSet(GATKSVVCFConstants.LOF,
            GATKSVVCFConstants.INT_EXON_DUP,
            GATKSVVCFConstants.DUP_PARTIAL,
            GATKSVVCFConstants.PARTIAL_EXON_DUP,
            GATKSVVCFConstants.COPY_GAIN,
            GATKSVVCFConstants.TSS_DUP);

    @VisibleForTesting
    protected enum StructuralVariantAnnotationType {
        DEL,
        DUP,
        INS,
        INV,
        CPX,
        BND,
        CTX,
        CNV
    }

    // Mini class to package SV type and interval into one object
    @VisibleForTesting
    protected static final class SVSegment {
        private final StructuralVariantAnnotationType intervalSVType;
        private final SimpleInterval interval;
        protected SVSegment(final StructuralVariantAnnotationType svType, final SimpleInterval interval) {
            this.intervalSVType = svType;
            this.interval = interval;
        }
        public StructuralVariantAnnotationType getIntervalSVType() {
            return intervalSVType;
        }
        public SimpleInterval getInterval() {
            return interval;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            final SVSegment svSegment = (SVSegment) o;
            return intervalSVType == svSegment.intervalSVType && interval.equals(svSegment.interval);
        }

        @Override
        public int hashCode() {
            return Objects.hash(intervalSVType, interval);
        }
    }

    // Container class for all SVIntervalTree trees created from the GTF
    @VisibleForTesting
    public static final class GTFIntervalTreesContainer {
        private final SVIntervalTree<GencodeGtfTranscriptFeature> transcriptIntervalTree;
        private final SVIntervalTree<String> promoterIntervalTree;
        private final SVIntervalTree<String> transcriptionStartSiteTree;
        protected GTFIntervalTreesContainer(final SVIntervalTree<GencodeGtfTranscriptFeature> transcriptIntervalTree,
                                            final SVIntervalTree<String> promoterIntervalTree,
                                            final SVIntervalTree<String> transcriptionStartSiteTree) {
            this.transcriptIntervalTree = transcriptIntervalTree;
            this.promoterIntervalTree = promoterIntervalTree;
            this.transcriptionStartSiteTree = transcriptionStartSiteTree;
        }

        public SVIntervalTree<GencodeGtfTranscriptFeature> getTranscriptIntervalTree() {
            return transcriptIntervalTree;
        }

        public SVIntervalTree<String> getPromoterIntervalTree() {
            return promoterIntervalTree;
        }

        public SVIntervalTree<String> getTranscriptionStartSiteTree() {
            return transcriptionStartSiteTree;
        }
    }

    public SVAnnotateEngine(final GTFIntervalTreesContainer gtfIntervalTrees,
                            final SVIntervalTree<String> nonCodingIntervalTree,
                            final SAMSequenceDictionary sequenceDictionary,
                            final int maxBreakendLen) {
        this.gtfIntervalTrees = gtfIntervalTrees;
        this.nonCodingIntervalTree = nonCodingIntervalTree;
        this.sequenceDictionary = sequenceDictionary;
        this.maxBreakendLen = maxBreakendLen;
    }

    /**
     * Checks if a given variant overlaps the TSS of a given transcript
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - boolean True if variant overlaps TSS of transcript, False if not
     */
    private static boolean variantOverlapsTranscriptionStartSite(final SimpleInterval variantInterval,
                                                                 final GencodeGtfTranscriptFeature gtfTranscript) {
        final int tss = SVAnnotate.getTranscriptionStartSite(gtfTranscript);
        return variantInterval.overlaps(new SimpleInterval(gtfTranscript.getContig(), tss, tss));
    }

    /**
     * Checks if a given variant interval spans (contains) a given feature interval
     * @param variantInterval - SimpleInterval representing structural variant
     * @param featureInterval - SimpleInterval representing feature (ie. transcript, noncoding element)
     * @return - boolean True if variant spans feature, False otherwise
     */
    @VisibleForTesting
    protected static boolean variantSpansFeature(final SimpleInterval variantInterval,
                                                 final SimpleInterval featureInterval) {
        return variantInterval.contains(featureInterval);
    }

    /**
     * Counts variant endpoints (breakends) inside feature interval
     * @param variantInterval - SimpleInterval representing structural variant
     * @param featureInterval - SimpleInterval representing feature (ie. transcript, noncoding element)
     * @return - int count of variant breakends that fall within the feature interval (0, 1, or 2)
     */
    @VisibleForTesting
    protected static int countBreakendsInsideFeature(final SimpleInterval variantInterval,
                                                     final SimpleInterval featureInterval) {
        if (!featureInterval.overlaps(variantInterval) || variantInterval.contains(featureInterval)) {
            return 0;
        } else if (featureInterval.contains(variantInterval)) {
            return 2;
        } else {
            return 1;
        }
    }

    /**
     * Adds key:value pair to provided variant consequence dictionary
     * @param variantConsequenceDict - map of variant consequence -> feature name
     * @param key - key (consequence) to add
     * @param value - value (feature name) to add
     */
    @VisibleForTesting
    protected static void updateVariantConsequenceDict(final Map<String, Set<String>> variantConsequenceDict,
                                                       final String key, final String value) {
        variantConsequenceDict.putIfAbsent(key, new HashSet<>());
        variantConsequenceDict.get(key).add(value);
    }

    /**
     * Common functionality underlying protein-coding consequence annotation for insertion, deletion, breakend, or
     * inversion SVs
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - preliminary consequence of the variant on the transcript
     */
    private static String getSimpleConsequence(final SimpleInterval variantInterval,
                                               final GencodeGtfTranscriptFeature gtfTranscript) {
        final List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
        String consequence = GATKSVVCFConstants.INTRONIC;
        for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
            if (!variantInterval.overlaps(gtfFeature)) {
                continue;
            }
            if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.CDS) {
                return GATKSVVCFConstants.LOF;
            } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                consequence = GATKSVVCFConstants.UTR;
            }
        }
        return consequence;
    }


    /**
     * Get consequence of insertion variant on transcript
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - consequence of insertion variant on transcript
     */
    @VisibleForTesting
    protected static String annotateInsertion(final SimpleInterval variantInterval,
                                              final GencodeGtfTranscriptFeature gtfTranscript) {
        return getSimpleConsequence(variantInterval, gtfTranscript);
    }

    /**
     * Get consequence of deletion variant on transcript
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - consequence of deletion variant on transcript
     */
    @VisibleForTesting
    protected static String annotateDeletion(final SimpleInterval variantInterval,
                                             final GencodeGtfTranscriptFeature gtfTranscript) {
        // For DEL only, return LOF if the variant overlaps the transcription start site
        if (variantOverlapsTranscriptionStartSite(variantInterval, gtfTranscript)) {
            return GATKSVVCFConstants.LOF;
        } else {
            return getSimpleConsequence(variantInterval, gtfTranscript);
        }
    }

    /**
     * Get consequence of duplication variant on transcript
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - consequence of duplication variant on transcript
     */
    @VisibleForTesting
    protected static String annotateDuplication(final SimpleInterval variantInterval,
                                                final GencodeGtfTranscriptFeature gtfTranscript) {
        final SimpleInterval transcriptInterval = new SimpleInterval(gtfTranscript);
        if (variantSpansFeature(variantInterval, transcriptInterval)) {
            return GATKSVVCFConstants.COPY_GAIN;
        } else if (variantOverlapsTranscriptionStartSite(variantInterval, gtfTranscript)) {
            return GATKSVVCFConstants.TSS_DUP;
        } else if (!transcriptInterval.contains(variantInterval)) {
            return GATKSVVCFConstants.DUP_PARTIAL;  // if one breakpoint inside transcript and one past the end
        } else {
            // both breakpoints inside transcript
            final List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
            int numBreakpointsInCDS = 0;
            int numBreakpointsInUTR = 0;
            int numCDSSpanned = 0;
            int numUTRSpanned = 0;
            for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
                if (!variantInterval.overlaps(gtfFeature)) {
                    continue;
                }
                final SimpleInterval featureInterval = new SimpleInterval(gtfFeature);
                if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.CDS) {
                    if (variantSpansFeature(variantInterval, featureInterval)) {
                        numCDSSpanned++;
                    } else {
                        numBreakpointsInCDS = numBreakpointsInCDS + countBreakendsInsideFeature(variantInterval, featureInterval);
                    }
                } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                    if (variantSpansFeature(variantInterval, featureInterval)) {
                        numUTRSpanned++;
                    } else {
                        numBreakpointsInUTR = numBreakpointsInUTR + countBreakendsInsideFeature(variantInterval, featureInterval);
                    }
                }
            }
            if (numBreakpointsInCDS == 2 || (numBreakpointsInCDS == 1 && numBreakpointsInUTR == 1)) {
                return GATKSVVCFConstants.LOF;
            } else if (numBreakpointsInCDS == 1) {
                return GATKSVVCFConstants.PARTIAL_EXON_DUP;
            } else if (numCDSSpanned > 0) {
                return GATKSVVCFConstants.INT_EXON_DUP;  // no breakpoints in CDS
            } else if (numBreakpointsInUTR > 0 || numUTRSpanned > 0) {
                return GATKSVVCFConstants.UTR;
            }
        }
        return GATKSVVCFConstants.INTRONIC;
    }

    /**
     * Get consequence for SV of type CNV (multiallelic copy number variant).
     * The consequence may depend on the individual's copy number at the site, so the variant is annotated as if it were
     * a biallelic duplication and then certain consequence categories are reclassified to PREDICTED_MSV_EXON_OVERLAP
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @param MSVExonOverlapClassifications - consequence categories from annotating a duplication that should be
     *                                      reclassified to PREDICTED_MSV_EXON_OVERLAP for a CNV
     * @return - consequence of CNV on transcript
     */
    @VisibleForTesting
    protected static String annotateCopyNumberVariant(final SimpleInterval variantInterval,
                                                      final GencodeGtfTranscriptFeature gtfTranscript,
                                                      final Set<String> MSVExonOverlapClassifications) {
        final String consequence = annotateDuplication(variantInterval, gtfTranscript);
        if (MSVExonOverlapClassifications.contains(consequence)) {
            return GATKSVVCFConstants.MSV_EXON_OVERLAP;
        } else {
            return consequence;
        }
    }

    /**
     * Get consequence of inversion variant on transcript
     * Shares common logic (simple consequence + TSS overlap -> LOF) and adds check for spanning inversions
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - consequence of inversion variant on transcript
     */
    @VisibleForTesting
    protected static String annotateInversion(final SimpleInterval variantInterval,
                                              final GencodeGtfTranscriptFeature gtfTranscript) {
        final SimpleInterval transcriptInterval = new SimpleInterval(gtfTranscript);
        if (variantSpansFeature(variantInterval, transcriptInterval)) {
            return GATKSVVCFConstants.INV_SPAN;
        } else {
            return annotateDeletion(variantInterval, gtfTranscript);
        }
    }

    /**
     * Get consequence of translocation on transcript
     * Only called with transcripts that overlap the translocation variant, so consequence is automatically LOF
     * because any translocation that breaks a gene is predicted to cause loss of function of that gene.
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - consequence of translocation on transcript
     */
    @VisibleForTesting
    protected static String annotateTranslocation(final SimpleInterval variantInterval,
                                                  final GencodeGtfTranscriptFeature gtfTranscript) {
        return GATKSVVCFConstants.LOF;
    }

    /**
     * Get consequence of breakend on transcript
     * Shares common logic with deletions (simple consequence + TSS overlap check) but low-confidence BNDs should not
     * be annotated as LOF, so LOF consequences are changed to BREAKEND_EXONIC
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - consequence of breakend on transcript
     */
    @VisibleForTesting
    protected static String annotateBreakend(final SimpleInterval variantInterval,
                                             final GencodeGtfTranscriptFeature gtfTranscript) {
        final String consequence = getSimpleConsequence(variantInterval, gtfTranscript);
        if (consequence.equals(GATKSVVCFConstants.LOF)) {
            return GATKSVVCFConstants.BREAKEND_EXON;
        }
        return consequence;
    }

    /**
     * Add consequence of structural variant on an overlapping transcript to consequence dictionary for variant
     * @param variantInterval - SimpleInterval representing structural variant
     * @param svType - SV type
     * @param transcript - protein-coding GTF transcript
     * @param variantConsequenceDict - running map of consequence -> feature name for variant to update
     */
    @VisibleForTesting
    protected void annotateTranscript(final SimpleInterval variantInterval,
                                             final StructuralVariantAnnotationType svType,
                                             final GencodeGtfTranscriptFeature transcript,
                                             final Map<String, Set<String>> variantConsequenceDict) {
        final String consequence;
        switch (svType) {
            case DEL:
                consequence = annotateDeletion(variantInterval, transcript);
                break;
            case INS:
                consequence = annotateInsertion(variantInterval, transcript);
                break;
            case DUP:
                consequence = annotateDuplication(variantInterval, transcript);
                break;
            case CNV:
                consequence = annotateCopyNumberVariant(variantInterval,transcript, MSV_EXON_OVERLAP_CLASSIFICATIONS);
                break;
            case INV:
                consequence = annotateInversion(variantInterval, transcript);
                break;
            case CTX:
                consequence = annotateTranslocation(variantInterval, transcript);
                break;
            case BND:
                consequence = annotateBreakend(variantInterval, transcript);
                break;
            default:
                consequence = null;
                break;
        }

        if (consequence != null) {
            updateVariantConsequenceDict(variantConsequenceDict, consequence, transcript.getGeneName());
        }
    }

    /**
     * Add annotations for any overlapping promoters for a variant to its consequence dictionary
     * @param variantInterval - SimpleInterval representing structural variant
     * @param variantConsequenceDict - running map of consequence -> feature name for variant to update
     */
    private void annotatePromoterOverlaps(final SimpleInterval variantInterval,
                                                 final Map<String, Set<String>> variantConsequenceDict) {
        final Set<String> codingAnnotationGenes = new HashSet<>();
        variantConsequenceDict.values().forEach(codingAnnotationGenes::addAll);
        final Iterator<SVIntervalTree.Entry<String>> promotersForVariant =
                gtfIntervalTrees.getPromoterIntervalTree().overlappers(
                        SVUtils.locatableToSVInterval(variantInterval, sequenceDictionary)
                );
        for (final Iterator<SVIntervalTree.Entry<String>> it = promotersForVariant; it.hasNext(); ) {
            final SVIntervalTree.Entry<String> promoterEntry = it.next();
            final String promoterName = promoterEntry.getValue();
            if (!codingAnnotationGenes.contains(promoterName)) {
                updateVariantConsequenceDict(variantConsequenceDict, GATKSVVCFConstants.PROMOTER, promoterName);
            }
        }
    }

    /**
     * Add annotations for any overlapping noncoding elements for a variant to its consequence dictionary
     * @param variantInterval - SimpleInterval representing structural variant
     * @param variantConsequenceDict - running map of consequence -> feature name for variant to update
     */
    private void annotateNonCodingOverlaps(final SimpleInterval variantInterval,
                                                  final Map<String, Set<String>> variantConsequenceDict) {
        final Iterator<SVIntervalTree.Entry<String>> nonCodingFeaturesForVariant =
                nonCodingIntervalTree.overlappers(SVUtils.locatableToSVInterval(variantInterval, sequenceDictionary));
        for (Iterator<SVIntervalTree.Entry<String>> it = nonCodingFeaturesForVariant; it.hasNext(); ) {
            SVIntervalTree.Entry<String> featureEntry = it.next();
            final String consequence =
                    variantSpansFeature(variantInterval, featureEntry.getInterval().toSimpleInterval(sequenceDictionary)) ?
                            GATKSVVCFConstants.NONCODING_SPAN : GATKSVVCFConstants.NONCODING_BREAKPOINT;
            updateVariantConsequenceDict(variantConsequenceDict, consequence, featureEntry.getValue());
        }
    }

    /**
     * Add nearest TSS annotation to a variant's consequence dictionary
     * To be run on intergenic variants that don't overlap promoters
     * @param variantInterval - SimpleInterval representing structural variant
     * @param variantConsequenceDict - running map of consequence -> feature name for variant to update
     */
    @VisibleForTesting
    protected void annotateNearestTranscriptionStartSite(final SimpleInterval variantInterval,
                                                         final Map<String, Set<String>> variantConsequenceDict) {
        // TODO: keep all nearest TSS for dispersed CPX / CTX or choose closest?
        final int variantContigID = SVUtils.getContigIDFromName(variantInterval.getContig(), sequenceDictionary);
        final SVInterval svInterval = SVUtils.locatableToSVInterval(variantInterval, sequenceDictionary);
        final SVIntervalTree.Entry<String> nearestBefore = gtfIntervalTrees.getTranscriptionStartSiteTree().max(svInterval);
        final SVIntervalTree.Entry<String> nearestAfter = gtfIntervalTrees.getTranscriptionStartSiteTree().min(svInterval);
        // nearest TSS only "valid" for annotation if non-null and on the same contig as the variant
        final boolean beforeValid = nearestBefore != null && nearestBefore.getInterval().getContig() == variantContigID;
        final boolean afterValid = nearestAfter != null && nearestAfter.getInterval().getContig() == variantContigID;
        // only update if at least one TSS is valid
        if (beforeValid || afterValid) {
            // set distance to closest valid TSS
            final int distanceBefore = beforeValid ? nearestBefore.getInterval().gapLen(svInterval) : Integer.MAX_VALUE;
            final int distanceAfter = afterValid ? svInterval.gapLen(nearestAfter.getInterval()) : Integer.MAX_VALUE;
            final String nearestTSSGeneName =
                    (distanceBefore < distanceAfter) ? nearestBefore.getValue() : nearestAfter.getValue();
            updateVariantConsequenceDict(variantConsequenceDict, GATKSVVCFConstants.NEAREST_TSS, nearestTSSGeneName);
        }
    }

    /**
     * Get SV type for variant from ALT field of VCF
     * @param variant - VCF record
     * @return - SV type determined from ALT field of VCF
     */
    @VisibleForTesting
    protected static StructuralVariantAnnotationType getSVType(final VariantContext variant) {
        if (variant.getAlternateAlleles().size() > 1) {
            throw new IllegalArgumentException("Expected single ALT allele, found multiple: " +
                    variant.getAlternateAlleles());
        }
        final Allele alt = variant.getAlternateAllele(0);
        if (alt.isBreakpoint()) {
            if (variant.hasAttribute(GATKSVVCFConstants.CPX_INTERVALS)) {
                return StructuralVariantAnnotationType.CPX;
            }
            return StructuralVariantAnnotationType.BND;
        } else if (alt.isSymbolic()) {
            if (alt.toString().contains("INS")) {
                // account for <INS:ME>, etc. types
                return StructuralVariantAnnotationType.INS;
            } else {
                // parse ALT as symbolic allele, assuming format <SVTYPE>
                return StructuralVariantAnnotationType.valueOf(alt.toString().substring(1, alt.toString().length()-1));
            }
        } else {
            throw new IllegalArgumentException("Unexpected ALT allele: " + alt +
                    ". Expected breakpoint or symbolic ALT allele representing a structural variant record.");
        }
    }

    /**
     * Add protein-coding annotations for any transcripts overlapping the variant to the variant consequence dictionary
     * @param variantInterval - SimpleInterval representing structural variant
     * @param svType - SV type
     * @param variantConsequenceDict - running map of consequence -> feature name for variant to update
     */
    @VisibleForTesting
    protected void annotateGeneOverlaps(final SimpleInterval variantInterval,
                                               final StructuralVariantAnnotationType svType,
                                               final Map<String, Set<String>> variantConsequenceDict) {
        final Iterator<SVIntervalTree.Entry<GencodeGtfTranscriptFeature>> gtfTranscriptsForVariant =
                gtfIntervalTrees.getTranscriptIntervalTree().overlappers(
                        SVUtils.locatableToSVInterval(variantInterval, sequenceDictionary)
                );
        for (Iterator<SVIntervalTree.Entry<GencodeGtfTranscriptFeature>> it = gtfTranscriptsForVariant; it.hasNext(); ) {
            SVIntervalTree.Entry<GencodeGtfTranscriptFeature> transcriptEntry = it.next();
            annotateTranscript(variantInterval, svType, transcriptEntry.getValue(), variantConsequenceDict);
        }
    }

    /**
     * Parse one interval string from CPX_INTERVALS INFO field into an SVSegment representing the SV type and
     * interval of one of the components of the complex event
     * @param cpxInterval - one element from CPX_INTERVALS list, a string representing one component of complex SV
     * @return - SVSegment representing one component of the complex SV (type and interval)
     */
    @VisibleForTesting
    protected static SVSegment parseCPXIntervalString(final String cpxInterval) {
        final String[] parsed = cpxInterval.split("_");
        final StructuralVariantAnnotationType svTypeForInterval = StructuralVariantAnnotationType.valueOf(parsed[0]);
        final SimpleInterval interval = new SimpleInterval(parsed[1]);
        return new SVSegment(svTypeForInterval, interval);
    }

    /**
     * Get SV type to use for annotation for a breakend VCF record
     * Breakend may represent BND, CTX, or DEL / DUP if the user specifies {@code SVAnnotate.MAX_BND_LEN_NAME}
     * @param variant - SimpleInterval representing structural variant
     * @param complexType - type of complex event - CPX_TYPE INFO field value
     * @param maxBreakendLen - Max size of BND in bp to annotate as DEL / DUP if applicable
     * @param svLen - SV length in bp - SVLEN INFO field value
     * @param chrom - chromosome where the variant is located - CHROM field of VCF
     * @param chr2 - second chromosome for the variant - CHR2 INFO field value
     * @return - SV type to use for annotation of breakend record
     */
    private static StructuralVariantAnnotationType getAnnotationTypeForBreakend(final VariantContext variant,
                                                                                final String complexType,
                                                                                final int maxBreakendLen,
                                                                                final int svLen,
                                                                                final String chrom, final String chr2) {
        if (complexType != null && complexType.contains("CTX")) {
            return StructuralVariantAnnotationType.CTX;
        } else if (maxBreakendLen > 0 && chr2 != null && chrom.equals(chr2) && svLen <= maxBreakendLen) {
            // if maxBreakendLenForOverlapAnnotation argument provided, annotate as DUP or DEL if applicable
            final String strand = variant.getAttributeAsString(GATKSVVCFConstants.STRANDS_ATTRIBUTE, null);
            if (strand == null) {
                return StructuralVariantAnnotationType.BND;  // not enough info to annotate as DEL or DUP TODO: throw error?
            }
            if (strand.equals(GATKSVVCFConstants.BND_DELETION_STRANDS)) {
                return StructuralVariantAnnotationType.DEL;
            } else if (strand.equals(GATKSVVCFConstants.BND_DUPLICATION_STRANDS)) {
                return StructuralVariantAnnotationType.DUP;
            }
        }
        return StructuralVariantAnnotationType.BND;
    }

    /**
     * Get list of SVSegments (type and interval) to annotate as one variant for a given VCF record
     * For simple variants, the list will contain one element representing the SV. Exceptions include:
     * - reciprocal translocations yield two positions (one on each chromosome)
     * - complex SVs yield multiple segments, one for each component of the complex event
     * - breakend records containing two breakpoints will yield two positions
     * @param variant - VCF record
     * @param overallSVType - SV type determined from ALT field of VCF record
     * @param maxBreakendLen - Max size of BND in bp to annotate as DEL / DUP if applicable
     * @return - list of SVSegments (type and interval) to annotate as part of the variant
     */
    @VisibleForTesting
    protected static List<SVSegment> getSVSegments(final VariantContext variant,
                                                   final StructuralVariantAnnotationType overallSVType,
                                                   final int maxBreakendLen) {
        final List<SVSegment> intervals;
        final String complexType = variant.getAttributeAsString(GATKSVVCFConstants.CPX_TYPE, null);
        final String chrom = variant.getContig();
        final int pos = variant.getStart();
        final String chr2 = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
        final int end2 = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, pos);
        if (overallSVType.equals(StructuralVariantAnnotationType.CPX)) {
            final List<String> cpxIntervalsString = variant.getAttributeAsStringList(GATKSVVCFConstants.CPX_INTERVALS, null);
            if (cpxIntervalsString == null) {
                throw new UserException("Complex (CPX) variant must contain CPX_INTERVALS INFO field");
            }
            if (complexType == null) {
                throw new UserException("Complex (CPX) variant must contain CPX_TYPE INFO field");
            }
            intervals = new ArrayList<>(cpxIntervalsString.size() + 1);
            for (final String cpxInterval : cpxIntervalsString) {
                intervals.add(parseCPXIntervalString(cpxInterval));
            }
            if (complexType.contains("dDUP")) {
                intervals.add(new SVSegment(StructuralVariantAnnotationType.INS,
                        new SimpleInterval(chrom, pos, pos + 1)));
            }
        } else if (overallSVType.equals(StructuralVariantAnnotationType.CTX)) {
            intervals = new ArrayList<>(2);
            intervals.add(new SVSegment(overallSVType, new SimpleInterval(variant)));  // CHROM:POS-POS+1
            // annotate both breakpoints of translocation - CHR2:END2-END2+1
            if (chr2 == null) {
                throw new UserException("Translocation (CTX) variant represented as a single record must contain CHR2 INFO field");
            }
            intervals.add(new SVSegment(overallSVType,
                    new SimpleInterval(chr2, end2, end2 + 1)));
        } else if (overallSVType.equals(StructuralVariantAnnotationType.BND)){
            intervals = new ArrayList<>(2);
            final int svLen = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
            // if BND representation of CTX event, get intervals as if BND but annotate as CTX
            final StructuralVariantAnnotationType annotateAs = getAnnotationTypeForBreakend(variant, complexType,
                    maxBreakendLen, svLen, chrom, chr2);
            if (annotateAs.equals(StructuralVariantAnnotationType.DEL) ||
                    annotateAs.equals(StructuralVariantAnnotationType.DUP)) {
                intervals.add(new SVSegment(annotateAs, new SimpleInterval(chrom, pos, pos + svLen)));
            } else {
                intervals.add(new SVSegment(annotateAs, new SimpleInterval(chrom, pos, pos)));
                if (chr2 != null && chr2.equals(chrom)) {
                    // will either have END2 or SVLEN - check which is in INFO and use that for second breakpoint
                    if (svLen > 0) {
                        intervals.add(new SVSegment(annotateAs,
                                new SimpleInterval(chrom, pos + svLen, pos + svLen)));
                    } else if (end2 != pos) {
                        intervals.add(new SVSegment(annotateAs, new SimpleInterval(chrom, end2, end2)));
                    }
                } else if (chr2 != null) {
                    // BND has CHR2 field and it is not equal to the main contig
                    intervals.add(new SVSegment(annotateAs, new SimpleInterval(chr2, end2, end2))); // if no end2 would just add duplicate segment which is ok
                } // TODO: else parse ALT field? Or just annotate this position and annotate the rest in the other records?
            }
        } else if (overallSVType.equals(StructuralVariantAnnotationType.INS)) {
            intervals = Collections.singletonList(new SVSegment(overallSVType,
                    new SimpleInterval(chrom, pos, pos + 1)));
        } else {
            intervals = Collections.singletonList(new SVSegment(overallSVType, new SimpleInterval(variant)));
        }

        return intervals;
    }



    /**
     * Create a copy of the variant consequence dictionary in which the feature names for each consequence are sorted
     * in alphabetical order
     * @param variantConsequenceDict - running map of consequence -> feature name for variant to update
     * @return - a consequence -> feature name map in which the feature names are sorted alphabetically
     */
    @VisibleForTesting
    protected static Map<String, Object> sortVariantConsequenceDict(final Map<String,Set<String>> variantConsequenceDict) {
        final Map<String, Object> sorted = new HashMap<>();
        for (final String consequence : variantConsequenceDict.keySet()) {
            final List<String> sortedGenes = new ArrayList<>(variantConsequenceDict.get(consequence));
            Collections.sort(sortedGenes);
            sorted.put(consequence, sortedGenes);
        }
        return sorted;
    }

    /**
     * Create a consequence -> feature name map and add all annotations for protein-coding, promoter, nearest TSS,
     * and noncoding consequences for a variant
     * @param variant - VCF record
     * @return - map of consequence -> feature name containing all annotations for the variant
     */
    @VisibleForTesting
    protected Map<String, Object> annotateStructuralVariant(final VariantContext variant) {
        final Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        final StructuralVariantAnnotationType overallSVType = getSVType(variant);
        final List<SVSegment> svSegments = getSVSegments(variant, overallSVType, maxBreakendLen);
        if (gtfIntervalTrees != null && gtfIntervalTrees.getTranscriptIntervalTree() != null) {
            for (SVSegment svSegment : svSegments) {
                annotateGeneOverlaps(svSegment.getInterval(), svSegment.getIntervalSVType(), variantConsequenceDict);
            }
        }

        // if variant consequence dictionary is empty (no protein-coding annotations), apply INTERGENIC flag
        final boolean noCodingAnnotations = variantConsequenceDict.isEmpty();

        // then annotate promoter overlaps and non-coding feature overlaps
        if (gtfIntervalTrees != null && gtfIntervalTrees.getPromoterIntervalTree() != null) {
            for (final SVSegment svSegment : svSegments) {
                annotatePromoterOverlaps(svSegment.getInterval(), variantConsequenceDict);
            }
        }

        if (nonCodingIntervalTree != null) {
            for (SVSegment svSegment : svSegments) {
                annotateNonCodingOverlaps(svSegment.getInterval(), variantConsequenceDict);
            }
        }

        // annotate nearest TSS for intergenic variants with no promoter overlaps
        if (gtfIntervalTrees != null && gtfIntervalTrees.getTranscriptionStartSiteTree() != null &&
                !variantConsequenceDict.containsKey(GATKSVVCFConstants.PROMOTER) && noCodingAnnotations) {
            for (SVSegment svSegment : svSegments) {
                annotateNearestTranscriptionStartSite(svSegment.getInterval(), variantConsequenceDict);
            }
        }

        final Map<String, Object> attributes = sortVariantConsequenceDict(variantConsequenceDict);
        if (gtfIntervalTrees != null && gtfIntervalTrees.getTranscriptIntervalTree() != null) {
            attributes.put(GATKSVVCFConstants.INTERGENIC, noCodingAnnotations);
        }
        return attributes;
    }

    /**
     * Create VariantContext for input variant with added functional annotation INFO keys
     * @param variant - input VCF record
     * @return - VariantContext equal to input + functional annotation INFO keys
     */
    public VariantContext createAnnotatedStructuralVariantContext(final VariantContext variant) {
        final Map<String, Object> attributes = annotateStructuralVariant(variant);
        return new VariantContextBuilder(variant)
                .putAttributes(attributes)
                .make();
    }


}
