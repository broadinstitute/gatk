package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
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
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.*;

import java.io.File;
import java.util.*;

/**
 * Adds gene overlap, predicted functional consequence, and noncoding element overlap annotations to
 * a structural variant (SV) VCF from the GATK-SV pipeline.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         VCF containing structural variant (SV) records from the GATK-SV pipeline
 *     </li>
 *     <li>
 *         GTF file containing primary or canonical transcripts (optional; required for protein-coding gene,
 *         transcription start site, and promoter overlap annotations)
 *     </li>
 *     <li>
 *         BED file of noncoding elements in which the fourth column specifies the type of element
 *         (optional; required for noncoding element overlap annotations)
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Annotated VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVAnnotate \
 *       -V structural.vcf.gz \
 *       --protein-coding-gtf canonical.gtf \
 *       --non-coding-bed noncoding.bed \
 *       -O annotated.vcf.gz
 * </pre>
 *
 * <h3>Annotation categories</h3>
 * <p>
 *     If a variant overlaps a gene, promoter, or noncoding element, the predicted functional impact will be annotated
 *     and the gene or noncoding element name listed. The list below describes the functional consequence categories
 *     in more detail. For complex variants (CPX), the annotations represent the union of the independent annotations
 *     of each component of the complex event according to its coordinates and SV type.
 * </p>
 *
 * <ul>
 *     <li><p><i>PREDICTED_LOF</i><br />
 *     Gene(s) on which the SV is predicted to have a loss-of-function (LOF) effect.
 *     The conditions for a LOF consequence depend on the SV type:
 *     <ul>
 *         <li>
 *             Deletion (DEL): the deletion overlaps any coding sequence (CDS) or the TSS
 *         </li>
 *         <li>
 *             Duplication (DUP): the duplication has both breakpoints in CDS, or one breakpoint in CDS and
 *             another in the 3' or 5' untranslated region (UTR)
 *         </li>
 *         <li>
 *             Insertion (INS): the insertion falls within CDS
 *         </li>
 *         <li>
 *             Inversion (INV): the inversion overlaps any coding sequence (CDS) or the TSS, except if it spans
 *             the entire gene (<i>PREDICTED_INV_SPAN</i>)
 *         </li>
 *         <li>
 *             Translocation (CTX): any translocation breakpoint falls within the transcript
 *         </li>
 *         <li>
 *             Multiallelic copy number variant (CNV): not annotated as <i>PREDICTED_LOF</i>.
 *             See <i>PREDICTED_MSV_EXON_OVERLAP</i>
 *         </li>
 *         <li>
 *             Breakend (BND): not annotated as <i>PREDICTED_LOF</i>. See <i>PREDICTED_BREAKEND_EXONIC</i>
 *         </li>
 *     </ul>
 *     </p></li>
 *     <li><p><i>PREDICTED_COPY_GAIN</i><br />
 *     Gene(s) on which the SV is predicted to have a copy-gain effect. This occurs when a duplication spans the entire
 *     transcript, from the first base of the 5' UTR to the last base of the 3' UTR. </p></li>
 *     <li><p><i>PREDICTED_INTRAGENIC_EXON_DUP</i><br />
 *     Gene(s) on which the SV is predicted to result in intragenic exonic duplication without breaking any coding
 *     sequences. This occurs when a duplication spans at least one coding exon and neither breakpoint is in CDS; both
 *     breakpoints are in UTR or intron. The result is that intact exons are duplicated within the boundaries of the
 *     gene body. </p></li>
 *     <li><p><i>PREDICTED_PARTIAL_EXON_DUP</i><br />
 *     Gene(s) where the duplication SV has one breakpoint in the coding sequence. This occurs when a duplication
 *     has exactly one breakpoint in CDS and the other breakpoint is in intron or UTR. When the duplication is in
 *     tandem, the result is that the endogenous copy of the breakpoint-harboring exon remains intact and a partial
 *     duplicate copy of that exon is also found elsewhere in the same gene.</p></li>
 *     <li><p><i>PREDICTED_TSS_DUP</i><br />
 *     Gene(s) for which the SV is predicted to duplicate the transcription start site (TSS). This occurs when a
 *     duplication has one breakpoint before the start of a transcript and the other breakpoint within the transcript.
 *     When the duplication is in tandem, the result is that there is one intact copy of the full endogenous gene, but
 *     an additional transcription start site is duplicated upstream (5’) of the endogenous TSS. </p></li>
 *     <li><p><i>PREDICTED_DUP_PARTIAL</i><br />
 *     Gene(s) which are partially overlapped by an SV's duplication, but the transcription start site is not
 *     duplicated. The partial duplication occurs when a duplication has one breakpoint within the transcript and one
 *     breakpoint after the end of the transcript. When the duplication is in tandem, the result is that there is one
 *     intact copy of the full endogenous gene.</p></li>
 *     <li><p><i>PREDICTED_INV_SPAN</i><br />
 *     Gene(s) which are entirely spanned by an SV's inversion. A whole-gene inversion occurs when an inversion spans
 *     the entire transcript, from the first base of the 5' UTR to the last base of the 3' UTR. </p></li>
 *     <li><p><i>PREDICTED_MSV_EXON_OVERLAP</i><br />
 *     Gene(s) on which the multiallelic CNV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN,
 *     DUP_PARTIAL, TSS_DUP, or PARTIAL_EXON_DUP annotation if the SV were biallelic. The functional impact of the
 *     multiallelic CNV on an individual sample depends on the copy number of the individual. </p></li>
 *     <li><p><i>PREDICTED_UTR</i><br />
 *     Gene(s) for which the SV is predicted to disrupt a UTR. This occurs when the SV has at least one breakpoint in
 *     a gene's 5' or 3' UTR but does not meet any of the criteria for a different gene-disrupting categorization
 *     above.</p></li>
 *     <li><p><i>PREDICTED_INTRONIC</i><br />
 *     Gene(s) where the SV was found to lie entirely within an intron. </p></li>
 *     <li><p><i>PREDICTED_BREAKEND_EXONIC</i><br />
 *     Gene(s) for which the SV breakend is predicted to fall in an exon. This category is reserved for breakend (BND)
 *     SVs with a breakpoint in CDS. </p></li>
 *     <li><p><i>PREDICTED_INTERGENIC</i><br />
 *     SV does not overlap any protein-coding gene transcripts in the GTF. </p></li>
 *     <li><p><i>PREDICTED_PROMOTER</i><br />
 *     Gene(s) for which the SV is predicted to overlap the promoter region. This occurs when the variant overlaps the
 *     predicted promoter region but does not overlap the transcript. The promoter region is inferred from the GTF as a
 *     window upstream of the TSS. The size of the window can be altered with the --promoter-window-length
 *     argument. </p></li>
 *     <li><p><i>PREDICTED_NEAREST_TSS</i><br />
 *     Nearest transcription start site to an intergenic variant. The gene with the nearest TSS to either side of the SV
 *     is annotated for intergenic variants that do not overlap any promoter regions. </p></li>
 *     <li><p><i>PREDICTED_NONCODING_SPAN</i><br />
 *     Class(es) of noncoding elements spanned by SV. </p></li>
 *     <li><p><i>PREDICTED_NONCODING_BREAKPOINT</i><br />
 *     Class(es) of noncoding elements disrupted by SV breakpoint. </p></li>
 *
 * </ul>
 */
@CommandLineProgramProperties(
        summary = "Adds predicted functional consequence, gene overlap, and noncoding element overlap annotations " +
                "to SV VCF from GATK-SV pipeline. " +
                "Input files are an SV VCF, a GTF file containing primary or canonical transcripts, " +
                "and a BED file containing noncoding elements. " +
                "Output file is an annotated SV VCF.",
        oneLineSummary = "Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class SVAnnotate extends VariantWalker {
    public static final String PROTEIN_CODING_GTF_NAME = "protein-coding-gtf";
    public static final String PROMOTER_WINDOW_NAME = "promoter-window-length";
    public static final String NON_CODING_BED_NAME = "non-coding-bed";
    public static final String MAX_BND_LEN_NAME = "max-breakend-as-cnv-length";

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file (if not provided, defaults to STDOUT)",
            common = false,
            optional = true
    )
    private GATKPath outputFile = null;

    @Argument(
            fullName=PROTEIN_CODING_GTF_NAME,
            doc="Protein-coding GTF file containing primary or canonical transcripts (1-2 transcripts per gene only)",
            optional=true
    )
    private File proteinCodingGTFFile;

    @Argument(
            fullName=PROMOTER_WINDOW_NAME,
            doc="Promoter window (bp) upstream of TSS. Promoters will be inferred as the {window} bases upstream of the TSS. Default: 1000",
            minValue=0, optional=true
    )
    private int promoterWindow = 1000;

    @Argument(
            fullName=NON_CODING_BED_NAME,
            doc="BED file (with header) containing non-coding features. Columns: chrom, start, end, name, score (.), strand",
            optional=true
    )
    private File nonCodingBedFile;

    @Argument(
            fullName=MAX_BND_LEN_NAME,
            doc="Length in bp. Provide to annotate BNDs smaller than this size as deletions or duplications if applicable. Recommended value: < 2000000",
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
    private SAMSequenceDictionary sequenceDictionary;

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
    protected static final class GTFIntervalTreesContainer {
        protected final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree;
        protected final SVIntervalTree<String> promoterIntervalTree;
        protected final SVIntervalTree<String> transcriptionStartSiteTree;
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
        // get contigs from VCF
        sequenceDictionary = header.getSequenceDictionary();
        // TODO: more elegant way to make reference inputs optional than checking 10x?
        // Load protein-coding GTF data into memory as interval tree of transcripts if GTF provided
        if (proteinCodingGTFFile != null) {
            final FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource = new FeatureDataSource<>(proteinCodingGTFFile);
            final GTFIntervalTreesContainer gtfTrees = buildIntervalTreesFromGTF(proteinCodingGTFSource, sequenceDictionary, promoterWindow);
            gtfIntervalTree = gtfTrees.gtfIntervalTree;
            promoterIntervalTree = gtfTrees.promoterIntervalTree;
            transcriptionStartSiteTree = gtfTrees.transcriptionStartSiteTree;
        }

        // Load noncoding BED file into memory as interval tree of noncoding elements if BED provided
        if (nonCodingBedFile != null) {
            final FeatureDataSource<FullBEDFeature> nonCodingSource = new FeatureDataSource<>(nonCodingBedFile);
            nonCodingIntervalTree = buildIntervalTreeFromBED(nonCodingSource, sequenceDictionary);
        }

        vcfWriter = createVCFWriter(outputFile);
        updateAndWriteHeader(header);
    }

    /**
     * Check if strand of transcript from GTF is negative
     * @param transcript - protein-coding GTF transcript to check
     * @return - boolean True if negative strand, False if positive
     */
    @VisibleForTesting
    protected static boolean isNegativeStrand(final GencodeGtfTranscriptFeature transcript) {
        return transcript.getGenomicStrand().equals(Strand.decode("-"));
    }

    /**
     * Get transcription start site from GTF transcript: first (5') base of transcript
     * @param transcript - protein-coding GTF transcript
     * @return - int position of TSS on transcript (start of transcript if + strand, end if - strand)
     */
    @VisibleForTesting
    protected static int getTranscriptionStartSite(final GencodeGtfTranscriptFeature transcript) {
        final boolean negativeStrand = isNegativeStrand(transcript);
        return negativeStrand ? transcript.getEnd() : transcript.getStart();  // GTF codec reassigns start to be earlier in contig;
    }

    /**
     * Get promoter interval for GTF transcript: user-defined window upstream (5') of TSS
     * @param transcript - protein-coding GTF transcript
     * @param promoterWindow - size of promoter window in bp
     * @return - SimpleInterval representing promoter region of size promoterWindow upstream of TSS
     */
    @VisibleForTesting
    protected static SimpleInterval getPromoterInterval(final GencodeGtfTranscriptFeature transcript,
                                                        final int promoterWindow) {
        final int tss = getTranscriptionStartSite(transcript);
        final boolean negativeStrand = isNegativeStrand(transcript);
        final int promoterLeft = negativeStrand ? tss + 1 : Math.max(tss - promoterWindow, 1);
        final int promoterRight = negativeStrand ? tss + promoterWindow : tss - 1;
        return new SimpleInterval(transcript.getContig(), promoterLeft, promoterRight);
    }

    /**
     * Checks if a given variant overlaps the TSS of a given transcript
     * @param variantInterval - SimpleInterval representing structural variant
     * @param gtfTranscript - protein-coding GTF transcript
     * @return - boolean True if variant overlaps TSS of transcript, False if not
     */
    private static boolean variantOverlapsTranscriptionStartSite(final SimpleInterval variantInterval,
                                                                 final GencodeGtfTranscriptFeature gtfTranscript) {
        final int tss = getTranscriptionStartSite(gtfTranscript);
        return variantInterval.overlaps(new SimpleInterval(gtfTranscript.getContig(), tss, tss));
    }

    /**
     * Builds transcript, TSS, and promoter interval trees from protein-coding GTF
     * @param proteinCodingGTFSource - GTF as FeatureDataSource
     * @param sequenceDictionary - SAMSequenceDictionary for VCF
     * @param promoterWindow - size of promoter window in bp
     * @return - container class packaging transcript, promoter, and TSS interval trees for annotation
     */
    @VisibleForTesting
    protected static GTFIntervalTreesContainer buildIntervalTreesFromGTF(final FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource,
                                                                final SAMSequenceDictionary sequenceDictionary, final int promoterWindow) {
        final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree = new SVIntervalTree<>();
        final SVIntervalTree<String> promoterIntervalTree = new SVIntervalTree<>();
        final SVIntervalTree<String> transcriptionStartSiteTree = new SVIntervalTree<>();
        for (final GencodeGtfGeneFeature gene : proteinCodingGTFSource) {
            final List<GencodeGtfTranscriptFeature> transcriptsForGene = gene.getTranscripts();
            for (GencodeGtfTranscriptFeature transcript : transcriptsForGene) {
                final int contigID = sequenceDictionary.getSequenceIndex(transcript.getContig());
                if (contigID < 0) {
                    continue; // if GTF input contains chromosome not in VCF sequence dictionary, just ignore it
                }
                gtfIntervalTree.put(SVUtils.locatableToSVInterval(transcript, sequenceDictionary), transcript);
                final String geneName = transcript.getGeneName();
                final int tss = getTranscriptionStartSite(transcript);
                transcriptionStartSiteTree.put(new SVInterval(contigID, tss, tss + 1), geneName);
                SimpleInterval promoterInterval = getPromoterInterval(transcript, promoterWindow);
                promoterIntervalTree.put(SVUtils.locatableToSVInterval(promoterInterval, sequenceDictionary), geneName);
            }
        }
        return new GTFIntervalTreesContainer(gtfIntervalTree, promoterIntervalTree, transcriptionStartSiteTree);
    }

    /**
     * Builds interval tree of noncoding elements to annotate from BED file input
     * @param BEDSource - noncoding element BED file as FeatureDataSource
     * @param sequenceDictionary - SAMSequenceDictionary for VCF
     * @return - SVIntervalTree of nonocoding elements for annotation
     */
    @VisibleForTesting
    protected static SVIntervalTree<String> buildIntervalTreeFromBED(final FeatureDataSource<FullBEDFeature> BEDSource,
                                                            final SAMSequenceDictionary sequenceDictionary) {
        final SVIntervalTree<String> BEDIntervalTree = new SVIntervalTree<>();
        for (final FullBEDFeature feature : BEDSource) {
            // BED feature class already does start+1 conversion to 1-based closed interval
            try {
                BEDIntervalTree.put(SVUtils.locatableToSVInterval(feature, sequenceDictionary), feature.getName());
            } catch (IllegalArgumentException e) {
                continue;  // if BED input contains chromosome not in VCF sequence dictionary, just ignore it
            }
        }
        return BEDIntervalTree;
    }

    /**
     * Adds SV functional annotation INFO keys to VCF header
     * @param header - starting VCF header to which to add INFO keys
     */
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
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MSV_EXON_OVERLAP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the multiallelic SV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN, DUP_PARTIAL, TSS_DUP, or PARTIAL_EXON_DUP annotation if the SV were biallelic."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PROMOTER, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to overlap the promoter region."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BREAKEND_EXON, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV breakend is predicted to fall in an exon."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INTERGENIC, 0, VCFHeaderLineType.Flag, "SV does not overlap any protein-coding genes."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NONCODING_SPAN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements spanned by SV."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NONCODING_BREAKPOINT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements disrupted by SV breakpoint."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NEAREST_TSS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Nearest transcription start site to an intergenic variant."));

    }

    /**
     * Updates VCF header with SV annotation INFO keys and writes header
     * @param header - starting VCF header
     */
    private void updateAndWriteHeader(final VCFHeader header) {
        addAnnotationInfoKeysToHeader(header);
        vcfWriter.writeHeader(header);
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
     * @param MSVExonOverlapClassifications - consequence classes to reclassify for CNVs
     */
    @VisibleForTesting
    protected static void annotateTranscript(final SimpleInterval variantInterval,
                                             final StructuralVariantAnnotationType svType,
                                             final GencodeGtfTranscriptFeature transcript,
                                             final Map<String, Set<String>> variantConsequenceDict,
                                             final Set<String> MSVExonOverlapClassifications) {
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
     * @param promoterIntervalTree - interval tree of promoters from GTF to find promoters overlapping variant
     * @param sequenceDictionary - SAMSequenceDictionary from VCF
     */
    private static void annotatePromoterOverlaps(final SimpleInterval variantInterval,
                                                    final Map<String, Set<String>> variantConsequenceDict,
                                                    final SVIntervalTree<String> promoterIntervalTree,
                                                    final SAMSequenceDictionary sequenceDictionary) {
        final Set<String> codingAnnotationGenes = new HashSet<>();
        variantConsequenceDict.values().forEach(codingAnnotationGenes::addAll);
        final Iterator<SVIntervalTree.Entry<String>> promotersForVariant =
                promoterIntervalTree.overlappers(SVUtils.locatableToSVInterval(variantInterval, sequenceDictionary));
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
     * @param nonCodingIntervalTree - interval tree of noncoding elements from BED to find features overlapping variant
     * @param sequenceDictionary - SAMSequenceDictionary from VCF
     */
    private static void annotateNonCodingOverlaps(final SimpleInterval variantInterval,
                                                  final Map<String, Set<String>> variantConsequenceDict,
                                                  final SVIntervalTree<String> nonCodingIntervalTree,
                                                  final SAMSequenceDictionary sequenceDictionary) {
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
     * @param transcriptionStartSiteTree - interval tree of TSS locations from GTF to find nearest TSS to variant
     * @param sequenceDictionary - SAMSequenceDictionary from VCF
     */
    @VisibleForTesting
    protected static void annotateNearestTranscriptionStartSite(final SimpleInterval variantInterval,
                                                                final Map<String, Set<String>> variantConsequenceDict,
                                                                final SVIntervalTree<String> transcriptionStartSiteTree,
                                                                final SAMSequenceDictionary sequenceDictionary) {
        // TODO: keep all nearest TSS for dispersed CPX / CTX or choose closest?
        final int variantContigID = SVUtils.getContigIDFromName(variantInterval.getContig(), sequenceDictionary);
        final SVInterval svInterval = SVUtils.locatableToSVInterval(variantInterval, sequenceDictionary);
        final SVIntervalTree.Entry<String> nearestBefore = transcriptionStartSiteTree.max(svInterval);
        final SVIntervalTree.Entry<String> nearestAfter = transcriptionStartSiteTree.min(svInterval);
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
     * @param MSVExonOverlapClassifications - consequence classes to reclassify for CNVs
     * @param sequenceDictionary - SAMSequenceDictionary from VCF
     * @param gtfIntervalTree - interval tree of protein-coding transcripts from GTF to find overlappers with variant
     */
    @VisibleForTesting
    protected static void annotateGeneOverlaps(final SimpleInterval variantInterval,
                                               final StructuralVariantAnnotationType svType,
                                               final Map<String, Set<String>> variantConsequenceDict,
                                               final Set<String> MSVExonOverlapClassifications,
                                               final SAMSequenceDictionary sequenceDictionary,
                                               final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree) {
        final Iterator<SVIntervalTree.Entry<GencodeGtfTranscriptFeature>> gtfTranscriptsForVariant =
                gtfIntervalTree.overlappers(SVUtils.locatableToSVInterval(variantInterval, sequenceDictionary));
        for (Iterator<SVIntervalTree.Entry<GencodeGtfTranscriptFeature>> it = gtfTranscriptsForVariant; it.hasNext(); ) {
            SVIntervalTree.Entry<GencodeGtfTranscriptFeature> transcriptEntry = it.next();
            annotateTranscript(variantInterval, svType, transcriptEntry.getValue(), variantConsequenceDict, MSVExonOverlapClassifications);
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
     * Breakend may represent BND, CTX, or DEL / DUP if the user specifies {@value MAX_BND_LEN_NAME}
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
     * @param gtfIntervalTree - interval tree of protein-coding transcripts from GTF to find overlappers with variant
     * @param promoterIntervalTree - interval tree of promoters from GTF to find promoters overlapping variant
     * @param transcriptionStartSiteTree - interval tree of TSS locations from GTF to find nearest TSS to variant
     * @param nonCodingIntervalTree - interval tree of noncoding elements from BED to find features overlapping variant
     * @param MSVExonOverlapClassifications - consequence classes to reclassify for CNVs
     * @param sequenceDictionary - SAMSequenceDictionary from VCF
     * @param maxBreakendLen - Max size of BND in bp to annotate as DEL / DUP if applicable
     * @return - map of consequence -> feature name containing all annotations for the variant
     */
    @VisibleForTesting
    protected static Map<String, Object> annotateStructuralVariant(final VariantContext variant,
                                                              final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree,
                                                              final SVIntervalTree<String> promoterIntervalTree,
                                                              final SVIntervalTree<String> transcriptionStartSiteTree,
                                                              final SVIntervalTree<String> nonCodingIntervalTree,
                                                              final Set<String> MSVExonOverlapClassifications,
                                                              final SAMSequenceDictionary sequenceDictionary,
                                                              final int maxBreakendLen) {
        final Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        final StructuralVariantAnnotationType overallSVType = getSVType(variant);
        final List<SVSegment> svSegments = getSVSegments(variant, overallSVType, maxBreakendLen);
        if (gtfIntervalTree != null) {
            for (SVSegment svSegment : svSegments) {
                annotateGeneOverlaps(svSegment.getInterval(), svSegment.getIntervalSVType(), variantConsequenceDict,
                        MSVExonOverlapClassifications, sequenceDictionary, gtfIntervalTree);
            }
        }

        // if variant consequence dictionary is empty (no protein-coding annotations), apply INTERGENIC flag
        final boolean noCodingAnnotations = variantConsequenceDict.isEmpty();

        // then annotate promoter overlaps and non-coding feature overlaps
        if (promoterIntervalTree != null) {
            for (final SVSegment svSegment : svSegments) {
                annotatePromoterOverlaps(svSegment.getInterval(), variantConsequenceDict, promoterIntervalTree,
                        sequenceDictionary);
            }
        }

        if (nonCodingIntervalTree != null) {
            for (SVSegment svSegment : svSegments) {
                annotateNonCodingOverlaps(svSegment.getInterval(), variantConsequenceDict, nonCodingIntervalTree,
                        sequenceDictionary);
            }
        }

        // annotate nearest TSS for intergenic variants with no promoter overlaps
        if (transcriptionStartSiteTree != null && !variantConsequenceDict.containsKey(GATKSVVCFConstants.PROMOTER) &&
                noCodingAnnotations) {
            for (SVSegment svSegment : svSegments) {
                annotateNearestTranscriptionStartSite(svSegment.getInterval(), variantConsequenceDict,
                        transcriptionStartSiteTree, sequenceDictionary);
            }
        }

        final Map<String, Object> attributes = sortVariantConsequenceDict(variantConsequenceDict);
        if (gtfIntervalTree != null) {
            attributes.put(GATKSVVCFConstants.INTERGENIC, noCodingAnnotations);
        }
        return attributes;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        final Map<String, Object> attributes = annotateStructuralVariant(variant, gtfIntervalTree, promoterIntervalTree,
                transcriptionStartSiteTree, nonCodingIntervalTree, MSVExonOverlapClassifications, sequenceDictionary,
                maxBreakendLen);
        final VariantContext annotated = new VariantContextBuilder(variant)
                .putAttributes(attributes)
                .make();
        vcfWriter.add(annotated);
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
