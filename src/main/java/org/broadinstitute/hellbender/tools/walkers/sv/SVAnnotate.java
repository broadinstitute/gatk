package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.FullBEDFeature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.compress.utils.Sets;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.ClosedSVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.*;

import java.io.File;
import java.util.*;

import static java.util.Objects.isNull;
import static org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata.buildContigIDToNameArray;
import static org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata.buildContigNameToIDMap;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils.locatableToClosedSVInterval;

/**
 * Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline
 * Input files are an SV VCF and a GTF file containing primary or canonical transcripts
 * Output file is an annotated SV VCF
 */
@CommandLineProgramProperties(
        summary = "Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline. " +
                "Input files are an SV VCF, a GTF file containing primary or canonical transcripts, " +
                "and a BED file containing noncoding elements. " +
                "Output file is an annotated SV VCF.",
        oneLineSummary = "Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class SVAnnotate extends VariantWalker {

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file (if not provided, defaults to STDOUT)",
            common = false,
            optional = true
    )
    private GATKPath outputFile = null;

    @Argument(
            fullName="proteinCodingGTF",
            doc="protein-coding GTF file (canonical only)",
            optional=true
    )
    private File proteinCodingGTFFile;

    @Argument(
            fullName="promoterWindowLength",
            doc="Promoter window (bp) upstream of TSS. Promoters will be inferred as the {window} bases upstream of the TSS. Default: 1000",
            minValue=0, optional=true
    )
    private int promoterWindow = 1000;

    @Argument(
            fullName="nonCodingBed",
            doc="BED file (with header) containing non-coding features. Columns: chrom, start, end, name, score (.), strand",
            optional=true
    )
    private File nonCodingBedFile;

    @Argument(
            fullName="maxBreakendLenForOverlapAnnotation",
            doc="Length in bp. Provide to annotate BNDs smaller than this value as deletions or duplications if applicable. Recommended value: < 2000000",
            minValue=0, optional=true
    )
    private int maxBreakendLen = -1;

    private VariantContextWriter vcfWriter = null;
    private SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree;
    private SVIntervalTree<String> promoterIntervalTree;
    private SVIntervalTree<String> nonCodingIntervalTree;
    private SVIntervalTree<String> transcriptionStartSiteTree;
    private final Set<String> MSVExonOverlapClassifications = Sets.newHashSet(GATKSVVCFConstants.LOF,
                                                                                GATKSVVCFConstants.INT_EXON_DUP,
                                                                                GATKSVVCFConstants.DUP_PARTIAL,
                                                                                GATKSVVCFConstants.PARTIAL_EXON_DUP,
                                                                                GATKSVVCFConstants.COPY_GAIN,
                                                                                GATKSVVCFConstants.TSS_DUP);
    private Map<String, Integer> contigNameToID;
    private String[] contigIDToName;
    private int maxContigLength;
    private enum StructuralVariantAnnotationType {
        DEL,
        DUP,
        INS,
        INV,
        CPX,
        BND,
        CTX,
        CNV
    }

    private static Integer getContigIDFromName(String contigName, Map<String, Integer> contigNameToID) {
        try {
            Integer contigID = contigNameToID.get(contigName);
            return contigID;
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("Contig " + contigName + " not in provided contig ID to name map");
        }
    }

    // mini class for SV intervals (type and segment) within CPX events
    protected static final class SVSegment {
        protected final StructuralVariantAnnotationType intervalSVType;
        protected final SimpleInterval interval;
        private SVSegment(final StructuralVariantAnnotationType svType, final SimpleInterval interval) {
            this.intervalSVType = svType;
            this.interval = interval;
        }
    }

    protected static final class GTFIntervalTreesContainer {
        protected SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree;
        protected SVIntervalTree<String> promoterIntervalTree;
        protected SVIntervalTree<String> transcriptionStartSiteTree;
        private GTFIntervalTreesContainer(final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree,
                                          final SVIntervalTree<String> promoterIntervalTree,
                                          final SVIntervalTree<String> transcriptionStartSiteTree) {
            this.gtfIntervalTree = gtfIntervalTree;
            this.promoterIntervalTree = promoterIntervalTree;
            this.transcriptionStartSiteTree = transcriptionStartSiteTree;
        }
    }

    @Override
    public void onTraversalStart() {
        final VCFHeader header = getHeaderForVariants();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        maxContigLength = sequenceDictionary.getSequences().stream().mapToInt(SAMSequenceRecord::getSequenceLength).max().getAsInt() + 1;
        contigNameToID = buildContigNameToIDMap(sequenceDictionary);
        contigIDToName = buildContigIDToNameArray(contigNameToID);
        // TODO: more elegant way to make reference inputs optional than checking 10x?
        if (!isNull(proteinCodingGTFFile)) {
            final FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource = new FeatureDataSource<>(proteinCodingGTFFile);
            GTFIntervalTreesContainer gtfTrees = buildIntervalTreesFromGTF(proteinCodingGTFSource, contigNameToID, promoterWindow);
            gtfIntervalTree = gtfTrees.gtfIntervalTree;
            promoterIntervalTree = gtfTrees.promoterIntervalTree;
            transcriptionStartSiteTree = gtfTrees.transcriptionStartSiteTree;
        }

        if (!isNull(nonCodingBedFile)) {
            final FeatureDataSource<FullBEDFeature> nonCodingSource = new FeatureDataSource<>(nonCodingBedFile);
            nonCodingIntervalTree = buildIntervalTreeFromBED(nonCodingSource);
        }

        vcfWriter = createVCFWriter(outputFile);
        updateAndWriteHeader(header);
    }

    protected static boolean isNegativeStrand(GencodeGtfTranscriptFeature transcript) {
        return transcript.getGenomicStrand().equals(Strand.decode("-"));
    }

    protected static int getTranscriptionStartSite(GencodeGtfTranscriptFeature transcript) {
        final boolean negativeStrand = isNegativeStrand(transcript);
        final int tss = negativeStrand ? transcript.getEnd() : transcript.getStart();  // GTF codec reassigns start to be earlier in contig
        return tss;
    }

    protected static ClosedSVInterval getPromoterInterval(GencodeGtfTranscriptFeature transcript, int promoterWindow,
                                                          int contigID) {
        final int tss = getTranscriptionStartSite(transcript);
        final boolean negativeStrand = isNegativeStrand(transcript);
        final int promoterLeft = negativeStrand ? tss + 1 : Math.max(tss - promoterWindow, 1);
        final int promoterRight = negativeStrand ? tss + promoterWindow : tss - 1;
        return new ClosedSVInterval(contigID, promoterLeft, promoterRight);
    }

    private static boolean variantOverlapsTranscriptionStartSite(SimpleInterval variantInterval,
                                                                 GencodeGtfTranscriptFeature gtfTranscript) {
        final int tss = getTranscriptionStartSite(gtfTranscript);
        return variantInterval.overlaps(new SimpleInterval(gtfTranscript.getContig(), tss, tss));
    }

    protected static GTFIntervalTreesContainer buildIntervalTreesFromGTF(FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource,
                                                                Map<String,Integer> contigNameToID, int promoterWindow) {
        // builds GTF, TSS, and promoter interval trees
        SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree = new SVIntervalTree<>();
        SVIntervalTree<String> promoterIntervalTree = new SVIntervalTree<>();
        SVIntervalTree<String> transcriptionStartSiteTree = new SVIntervalTree<>();
        for (final GencodeGtfGeneFeature gene : proteinCodingGTFSource) {
            final List<GencodeGtfTranscriptFeature> transcriptsForGene = gene.getTranscripts();
            for (GencodeGtfTranscriptFeature transcript : transcriptsForGene) {
                final int contigID;
                try {
                    contigID = getContigIDFromName(transcript.getContig(), contigNameToID);
                } catch (IllegalArgumentException e) {
                    continue; // if GTF input contains chromosome not in VCF sequence dictionary, just ignore it
                }
                gtfIntervalTree.put(locatableToClosedSVInterval(transcript, contigNameToID), transcript);
                final String geneName = transcript.getGeneName();
                final int tss = getTranscriptionStartSite(transcript);
                transcriptionStartSiteTree.put(new ClosedSVInterval(contigID, tss, tss), geneName);
                ClosedSVInterval promoterInterval = getPromoterInterval(transcript, promoterWindow, contigID);
                promoterIntervalTree.put(promoterInterval, geneName);
            }
        }
        return new GTFIntervalTreesContainer(gtfIntervalTree, promoterIntervalTree, transcriptionStartSiteTree);
    }

    private SVIntervalTree<String> buildIntervalTreeFromBED(FeatureDataSource<FullBEDFeature> BEDSource) {
        SVIntervalTree<String> BEDIntervalTree = new SVIntervalTree<>();
        for (final FullBEDFeature feature : BEDSource) {
            // BED feature already does start+1 conversion to closed interval
            try {
                BEDIntervalTree.put(locatableToClosedSVInterval(feature, contigNameToID), feature.getName()); // TODO: handle this in a function to convert Locatable to SVInterval
            } catch (IllegalArgumentException e) {
                continue;  // if BED input contains chromosome not in VCF sequence dictionary, just ignore it
            }
        }
        return BEDIntervalTree;
    }


    private void addAnnotationInfoKeysToHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.LOF, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to have a loss-of-function effect."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INT_EXON_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to result in intragenic exonic duplication without breaking any coding sequences."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.COPY_GAIN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to have a copy-gain effect."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TSS_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to duplicate the transcription start site."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_PARTIAL, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) which are partially overlapped by an SV's duplication, but the transcription start site is not duplicated."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INTRONIC, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) where the SV was found to lie entirely within an intron."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PARTIAL_EXON_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) where the duplication SV has one breakpoint in the coding sequence."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INV_SPAN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) which are entirely spanned by an SV's inversion."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.UTR, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to disrupt a UTR."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MSV_EXON_OVERLAP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the multiallelic SV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN, DUP_PARTIAL, or PARTIAL_EXON_DUP annotation if the SV were biallelic."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PROMOTER, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to overlap the promoter region."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BREAKEND_EXON, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV breakend is predicted to fall in an exon."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INTERGENIC, 0, VCFHeaderLineType.Flag, "SV does not overlap coding sequence."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NONCODING_SPAN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements spanned by SV."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NONCODING_BREAKPOINT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements disrupted by SV breakpoint."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NEAREST_TSS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Nearest transcription start site to intragenic variants."));

    }

    private void updateAndWriteHeader(VCFHeader header) {
        addAnnotationInfoKeysToHeader(header);
        vcfWriter.writeHeader(header);
    }

    protected static boolean variantSpansFeature(final SimpleInterval variantInterval,
                                                 final SimpleInterval featureInterval) {
        return variantInterval.contains(featureInterval);
    }

    protected static int countBreakendsInsideFeature(final SimpleInterval variantInterval,
                                                     final SimpleInterval featureInterval) {
        int count = 0;
        if (variantInterval.getContig().equals(featureInterval.getContig())) {
            if (variantInterval.getStart() >= featureInterval.getStart() &&
                    variantInterval.getStart() <= featureInterval.getEnd()) {
                count++;
            }
            if (variantInterval.getEnd() >= featureInterval.getStart() &&
                    variantInterval.getEnd() <= featureInterval.getEnd()) {
                count++;
            }
        }
        return count;
    }

    protected static void updateVariantConsequenceDict(final Map<String, Set<String>> variantConsequenceDict,
                                                     final String key,
                                                     final String value) {
        variantConsequenceDict.putIfAbsent(key, new HashSet<>());
        variantConsequenceDict.get(key).add(value);
    }

    protected static String annotateInsertion(final SimpleInterval variantInterval,
                                               final GencodeGtfTranscriptFeature gtfTranscript) {
        final List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
        String consequence = GATKSVVCFConstants.INTRONIC;
        for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
            if (!variantInterval.overlaps(gtfFeature)) {
                continue;
            }
            if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.CDS) {
                consequence = GATKSVVCFConstants.LOF;
                break;
            } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                consequence = GATKSVVCFConstants.UTR;
            }
        }
        return consequence;
    }

    protected static String annotateDeletion(final SimpleInterval variantInterval,
                                    final GencodeGtfTranscriptFeature gtfTranscript) {
        // For DEL only, return LOF if the variant overlaps the transcription start site
        if (variantOverlapsTranscriptionStartSite(variantInterval, gtfTranscript)) {
            return GATKSVVCFConstants.LOF;
        } else {
            return annotateInsertion(variantInterval, gtfTranscript);
        }
    }

    protected static String annotateDuplication(final SimpleInterval variantInterval,
                                       final GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = GATKSVVCFConstants.INTRONIC;
        final SimpleInterval transcriptInterval = new SimpleInterval(gtfTranscript);
        if (variantSpansFeature(variantInterval, transcriptInterval)) {
            consequence = GATKSVVCFConstants.COPY_GAIN;
        } else if (variantOverlapsTranscriptionStartSite(variantInterval, gtfTranscript)) {
            consequence = GATKSVVCFConstants.TSS_DUP;
        } else if (countBreakendsInsideFeature(variantInterval, transcriptInterval) == 1) {
            consequence = GATKSVVCFConstants.DUP_PARTIAL;
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
                consequence = GATKSVVCFConstants.LOF;
            } else if (numBreakpointsInCDS == 1) {
                consequence = GATKSVVCFConstants.PARTIAL_EXON_DUP;
            } else if (numCDSSpanned > 0) {
                consequence = GATKSVVCFConstants.INT_EXON_DUP;  // no breakpoints in CDS
            } else if (numBreakpointsInUTR > 0 || numUTRSpanned > 0) {
                consequence = GATKSVVCFConstants.UTR;
            }
        }
        return consequence;
    }

    protected static String annotateCopyNumberVariant(final SimpleInterval variantInterval,
                                                      final GencodeGtfTranscriptFeature gtfTranscript,
                                                      final Set<String> MSVExonOverlapClassifications) {
        String consequence = annotateDuplication(variantInterval, gtfTranscript);
        if (MSVExonOverlapClassifications.contains(consequence)) {
            return GATKSVVCFConstants.MSV_EXON_OVERLAP;
        } else {
            return consequence;  // TODO: MCNV classifications ???
        }
    }

    protected static String annotateInversion(final SimpleInterval variantInterval,
                                     final GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence;
        final SimpleInterval transcriptInterval = new SimpleInterval(gtfTranscript);
        if (variantSpansFeature(variantInterval, transcriptInterval)) {
            consequence = GATKSVVCFConstants.INV_SPAN;
        } else {
            consequence = annotateDeletion(variantInterval, gtfTranscript);
        }
        return consequence;
    }

    protected static String annotateTranslocation(final SimpleInterval variantInterval,
                                         final GencodeGtfTranscriptFeature gtfTranscript) {
        // already checked for transcript overlap, and if a translocation breakpoint falls inside a gene it's automatically LOF
        // TODO: eliminate unnecessary function or keep for aesthetics/future flexibility?
        return GATKSVVCFConstants.LOF;
    }

    protected static String annotateBreakend(final SimpleInterval variantInterval,
                                    final GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = annotateInsertion(variantInterval, gtfTranscript);
        if (consequence.equals(GATKSVVCFConstants.LOF)) {
            consequence = GATKSVVCFConstants.BREAKEND_EXON;
        }
        return consequence;
    }

    protected static void annotateTranscript(final SimpleInterval variantInterval,
                                    final StructuralVariantAnnotationType svType,
                                    final GencodeGtfTranscriptFeature transcript,
                                    final Map<String, Set<String>> variantConsequenceDict,
                                    final Set<String> MSVExonOverlapClassifications) {
        String consequence = null;
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
                consequence = annotateCopyNumberVariant(variantInterval,transcript, MSVExonOverlapClassifications);
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
                break;
        }

        if (consequence != null) {
            updateVariantConsequenceDict(variantConsequenceDict, consequence, transcript.getGeneName());
        }
    }

    /**
     Annotate promoter overlaps and return a boolean: true if there are any promoter overlaps, false if not
     */
    private boolean annotatePromoterOverlaps(final SimpleInterval variantInterval,
                                             final Map<String, Set<String>> variantConsequenceDict) {
        boolean anyPromoterOverlaps = false;
        final Set<String> codingAnnotationGenes = new HashSet<>();
        variantConsequenceDict.values().forEach(codingAnnotationGenes::addAll);
        final Iterator<SVIntervalTree.Entry<String>> promotersForVariant =
                promoterIntervalTree.overlappers(locatableToClosedSVInterval(variantInterval, contigNameToID));
        for (Iterator<SVIntervalTree.Entry<String>> it = promotersForVariant; it.hasNext(); ) {
            SVIntervalTree.Entry<String> promoterEntry = it.next();
            String promoterName = promoterEntry.getValue();
            if (!codingAnnotationGenes.contains(promoterName)) {
                updateVariantConsequenceDict(variantConsequenceDict, GATKSVVCFConstants.PROMOTER, promoterName);
                anyPromoterOverlaps = true;
            }
        }
        return anyPromoterOverlaps;
    }

    private void annotateNonCodingOverlaps(final SimpleInterval variantInterval,
                                           final Map<String, Set<String>> variantConsequenceDict) {
        final Iterator<SVIntervalTree.Entry<String>> nonCodingFeaturesForVariant =
                nonCodingIntervalTree.overlappers(locatableToClosedSVInterval(variantInterval, contigNameToID));
        for (Iterator<SVIntervalTree.Entry<String>> it = nonCodingFeaturesForVariant; it.hasNext(); ) {
            SVIntervalTree.Entry<String> featureEntry = it.next();
            String consequence = GATKSVVCFConstants.NONCODING_BREAKPOINT;
            if (variantSpansFeature(variantInterval, featureEntry.getInterval().toSimpleInterval(contigIDToName))) {
                consequence = GATKSVVCFConstants.NONCODING_SPAN;
            }
            updateVariantConsequenceDict(variantConsequenceDict, consequence, featureEntry.getValue());
        }
    }


    protected static void annotateNearestTranscriptionStartSite(final SimpleInterval variantInterval,
                                                                final Map<String, Set<String>> variantConsequenceDict,
                                                                final SVIntervalTree<String> transcriptionStartSiteTree,
                                                                final int maxContigLength,
                                                                final int variantContigID) {
        // TODO: keep all nearest TSS for dispersed CPX / CTX or choose closest?
        final ClosedSVInterval svInterval = new ClosedSVInterval(variantContigID, variantInterval.getStart(), variantInterval.getEnd());
        final SVIntervalTree.Entry<String> nearestBefore = transcriptionStartSiteTree.max(svInterval);
        final SVIntervalTree.Entry<String> nearestAfter = transcriptionStartSiteTree.min(svInterval);
        // nearest TSS only "valid" for annotation if non-null and on the same contig as the variant
        final boolean beforeInvalid = (isNull(nearestBefore) || nearestBefore.getInterval().getContig() != variantContigID );
        final boolean afterInvalid = (isNull(nearestAfter) || nearestAfter.getInterval().getContig() != variantContigID );
        // if at least one result is valid, keep one with shorter distance
        if (!(beforeInvalid && afterInvalid)) {
            // if result is invalid, set distance to longest contig length so that other TSS will be kept
            final int distanceBefore = beforeInvalid ? maxContigLength : variantInterval.getStart() - nearestBefore.getInterval().getEnd();
            final int distanceAfter = afterInvalid ? maxContigLength : nearestAfter.getInterval().getStart() - variantInterval.getEnd();
            final String nearestTSSGeneName = (distanceBefore < distanceAfter) ? nearestBefore.getValue() : nearestAfter.getValue();
            updateVariantConsequenceDict(variantConsequenceDict, GATKSVVCFConstants.NEAREST_TSS, nearestTSSGeneName);
        }
        // TODO: return consequence instead? easier for unit tests...
    }

    private StructuralVariantAnnotationType getSVType(final VariantContext variant) {
        // TODO: haha majorly clean this up
        // return variant.getStructuralVariantType().name();
        final Allele alt = variant.getAlternateAllele(0); // TODO: any chance of multiallelic alt field for SV?
        if (alt.isBreakpoint()) {
            if (variant.hasAttribute(GATKSVVCFConstants.CPX_INTERVALS)) {
                return StructuralVariantAnnotationType.CPX;
            }
            return StructuralVariantAnnotationType.BND;
        } else if (alt.isSingleBreakend()) {
            throw new IllegalArgumentException("what even is single breakend??: " + alt);
        } else if (alt.isSymbolic()) {
            if (alt.toString().contains("INS")) {
                // account for <INS:ME>, etc. types
                return StructuralVariantAnnotationType.INS;
            } else {
                // parse ALT as symbolic allele, assuming format <SVTYPE>
                return StructuralVariantAnnotationType.valueOf(alt.toString().substring(1, alt.toString().length()-1));
            }
        } else {
            throw new IllegalArgumentException("Unexpected ALT allele: " + alt);
        }
    }

    protected static void annotateGeneOverlaps(final SimpleInterval variantInterval,
                                      final StructuralVariantAnnotationType svType,
                                      final Map<String, Set<String>> variantConsequenceDict,
                                      final Set<String> MSVExonOverlapClassifications,
                                      final Map<String,Integer> contigNameToID,
                                      final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree) {
        // TODO: method to convert SimpleInterval > SVInterval? Or switch to always using one or the other - could I use SimpleInterval as base for ClosedIntervalTree?
        final Iterator<SVIntervalTree.Entry<GencodeGtfTranscriptFeature>> gtfTranscriptsForVariant =
                gtfIntervalTree.overlappers(locatableToClosedSVInterval(variantInterval, contigNameToID));
        for (Iterator<SVIntervalTree.Entry<GencodeGtfTranscriptFeature>> it = gtfTranscriptsForVariant; it.hasNext(); ) {
            SVIntervalTree.Entry<GencodeGtfTranscriptFeature> transcriptEntry = it.next();
            annotateTranscript(variantInterval, svType, transcriptEntry.getValue(), variantConsequenceDict, MSVExonOverlapClassifications);
        }
    }

    protected static SVSegment parseCPXIntervalString(String cpxInterval) {
        final String[] parsed = cpxInterval.split("_");
        final StructuralVariantAnnotationType svTypeForInterval = StructuralVariantAnnotationType.valueOf(parsed[0]);
        final SimpleInterval interval = new SimpleInterval(parsed[1]);
        return new SVSegment(svTypeForInterval, interval);
    }

    private static StructuralVariantAnnotationType getAnnotationTypeForBreakend(VariantContext variant,
                                                                                int maxBreakendLen, int svLen,
                                                                                String chrom, String chr2) {
        StructuralVariantAnnotationType annotateAs = StructuralVariantAnnotationType.BND;
        if (variant.getAttributeAsString(GATKSVVCFConstants.CPX_TYPE, "NONE").contains("CTX")) {
            annotateAs = StructuralVariantAnnotationType.CTX;
        } else if (maxBreakendLen > 0 && chrom.equals(chr2) && svLen <= maxBreakendLen) {
            final String strand = variant.getAttributeAsString(GATKSVVCFConstants.STRANDS_ATTRIBUTE, "NONE");
            if (strand.equals(GATKSVVCFConstants.BND_DELETION_STRANDS)) {
                annotateAs = StructuralVariantAnnotationType.DEL;
            } else if (strand.equals(GATKSVVCFConstants.BND_DUPLICATION_STRANDS)) {
                annotateAs = StructuralVariantAnnotationType.DUP;
            }
        }
        return annotateAs;
    }

    private static List<SVSegment> getSVSegments(VariantContext variant, StructuralVariantAnnotationType overallSVType,
                                          int maxBreakendLen) {
        final List<SVSegment> intervals = new ArrayList<>();
        final String chrom = variant.getContig();
        final int pos = variant.getStart();
        final String chr2 = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, "NONE");
        final int end2 = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, pos); // TODO: what if no end2 for CTX, BND with 2nd break on other contig?
        if (overallSVType.equals(StructuralVariantAnnotationType.CPX)) {
            final List<String> cpxIntervalsString = variant.getAttributeAsStringList(GATKSVVCFConstants.CPX_INTERVALS, "NONE");
            for (String cpxInterval : cpxIntervalsString) {
                intervals.add(parseCPXIntervalString(cpxInterval));
            }
        } else if (overallSVType.equals(StructuralVariantAnnotationType.CTX)) {
            intervals.add(new SVSegment(overallSVType, new SimpleInterval(variant)));  // CHROM:POS-POS+1
            // annotate both breakpoints of translocation - CHR2:END2-END2+1
            intervals.add(new SVSegment(overallSVType,
                    new SimpleInterval(chr2, end2, end2 + 1))); // TODO: when to add +1 to end? CTX +1, BND without?
        } else if (overallSVType.equals(StructuralVariantAnnotationType.BND)){
            final int svLen = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0); // TODO: what if no SVLEN and chr2 == chrom?
            // if BND representation of CTX event, get intervals as if BND but annotate as CTX
            StructuralVariantAnnotationType annotateAs = getAnnotationTypeForBreakend(variant, maxBreakendLen, svLen, chrom, chr2);
            intervals.add(new SVSegment(annotateAs, new SimpleInterval(chrom, pos, pos)));
            if (chr2.equals(chrom)) {
                intervals.add(new SVSegment(annotateAs, new SimpleInterval(chrom, pos + svLen, pos + svLen)));
            } else if (!chr2.equals("NONE")) {
                // BND has CHR2 field and it is not equal to the main contig
                intervals.add(new SVSegment(annotateAs, new SimpleInterval(chr2, end2, end2))); // if no end2 would just add duplicate segment which is ok
            } // TODO: else parse ALT field? Or just annotate this position and annotate the rest in the other records?
        } else {
            intervals.add(new SVSegment(overallSVType, new SimpleInterval(variant)));
        }

        return intervals;
    }

    protected static Map<String,String> formatVariantConsequenceDict(Map<String,Set<String>> variantConsequenceDict) {
        Map<String,String> formatted = new HashMap<>();
        for (String consequence : variantConsequenceDict.keySet()) {
            List<String> sortedGenes = new ArrayList<>(variantConsequenceDict.get(consequence));
            Collections.sort(sortedGenes);
            formatted.put(consequence, String.join(",", sortedGenes));
        }
        return formatted;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        final Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        final StructuralVariantAnnotationType overallSVType = getSVType(variant);
        final List<SVSegment> svSegments = getSVSegments(variant, overallSVType, maxBreakendLen);
        if (!isNull(proteinCodingGTFFile)) {
            for (SVSegment svSegment : svSegments) {
                annotateGeneOverlaps(svSegment.interval, svSegment.intervalSVType, variantConsequenceDict,
                        MSVExonOverlapClassifications, contigNameToID, gtfIntervalTree);
            }
        }

        // if variant consequence dictionary is empty (no protein-coding annotations), apply INTERGENIC flag
        boolean noCodingAnnotations = variantConsequenceDict.isEmpty();

        // then annotate promoter overlaps and non-coding feature overlaps
        boolean anyPromoterOverlaps = false;
        if (!isNull(proteinCodingGTFFile)) {
            for (SVSegment svSegment : svSegments) {
                anyPromoterOverlaps = anyPromoterOverlaps || annotatePromoterOverlaps(svSegment.interval, variantConsequenceDict);
            }
        }

        if (!isNull(nonCodingBedFile)) {
            for (SVSegment svSegment : svSegments) {
                annotateNonCodingOverlaps(svSegment.interval, variantConsequenceDict);
            }
        }

        // annotate nearest TSS for intergenic variants with no promoter overlaps
        if (!isNull(proteinCodingGTFFile) && !anyPromoterOverlaps && noCodingAnnotations) {
            for (SVSegment svSegment : svSegments) {
                annotateNearestTranscriptionStartSite(svSegment.interval, variantConsequenceDict, transcriptionStartSiteTree,
                                                        maxContigLength, contigNameToID.get(svSegment.interval.getContig()));
            }
        }

        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        vcb.putAttributes(formatVariantConsequenceDict(variantConsequenceDict));
        if (!isNull(proteinCodingGTFFile)) {
            vcb.attribute(GATKSVVCFConstants.INTERGENIC, noCodingAnnotations);
        }
        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
