package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * A simple data object to hold a comparison between a reference sequence and an alternate allele.
 *
 * Created by jonn on 8/30/17.
 */
public class SequenceComparison {

    /**
     * Bases covering the region around the variant.
     * This includes bases that would overlap any of the variant or reference bases.
     * (i.e. if the reference allele is 'A' and the variant is 'GCGCG', this would be
     * 'ANNNNN' where 'N' is the correct base at that position in the reference sequence).
     * Contains {@link SequenceComparison#referenceWindow} bases of padding before and
     * after the total length.
     */
    private String referenceBases                     = null;

    /**
     * The number of bases in {@link SequenceComparison#referenceBases} before the
     * start of the reference Allele / variant.
     */
    private Integer referenceWindow                   = null;

    /**
     * The reference coding sequence for a the transcript of this sequence comparison.
     * This does NOT include introns.
     * Stored in the forward reading direction.  For NEGATIVE strand reads, must
     * reverse complement any bases retrieved.
     */
    private ReferenceSequence transcriptCodingSequence = null;

    /**
     * The reference coding sequence for a the transcript of this sequence comparison.
     * This does NOT include introns.
     * Stored in the forward reading direction.  For NEGATIVE strand reads, must
     * reverse complement any bases retrieved.
     */
    private ReferenceSequence referenceCodingSequence = null;

    /**
     * The contig on which this sequence comparison occurs.
     */
    private String  contig                               = null;

    /**
     * The strand on which this sequence comparison occurs.
     */
    private Strand strand                               = null;

    /**
     * The position (1-based, inclusive in genome coordinates relative to the start of
     * {@link SequenceComparison#contig}) of the start of the reference allele / variant.
     */
    private Integer alleleStart                          = null;

    /**
     * The position (1-based, inclusive) in transcript coordinates relative to the start of
     * the transcript of start of the reference allele / variant.
     */
    private Integer transcriptAlleleStart                = null;

    /**
     * The position (1-based, inclusive) in coding sequence coordinates relative to the start of
     * the coding region of the transcript of the start of the reference allele / variant.
     * This location is obtained by concatenating the exons together and counting from the start of
     * that sequence to where the reference allele / variant begins.
     */
    private Integer codingSequenceAlleleStart            = null;

    /**
     * The in-frame position (1-based, inclusive) in coding sequence coordinates relative to the start of
     * the coding region of the transcript of the start of the first codon containing the reference allele / variant.
     * This location is obtained by concatenating the exons together and counting from the start of
     * that sequence to where the reference allele / variant begins, then moving backwards to an in-frame
     * position, if necessary.
     */
    private Integer alignedCodingSequenceAlleleStart     = null;

    /**
     * The position (1-based, inclusive in genome coordinates - relative to the start of
     * {@link SequenceComparison#contig}) of the start of the exon that contains the
     * variant in this {@link SequenceComparison}.
     */
    private Integer exonStartPosition                    = null;

    /**
     * The position (1-based, inclusive in genome coordinates - relative to the start of
     * {@link SequenceComparison#contig}) of the end of the exon that contains the
     * variant in this {@link SequenceComparison}.
     */
    private Integer exonEndPosition                      = null;

    /**
     * The {@link String} representation of the reference allele, correctly represented for the {@link Strand} on which it appears.
     */
    private String referenceAllele                      = null;
    /**
     * An in-frame sequence of bases that overlaps the given reference allele based on the raw reference genome.
     * (i.e. This includes INTRONS.)
     */
    private String  alignedReferenceAllele               = null;
    /**
     * An in-frame sequence of bases that overlaps the given reference allele for the coding region only.
     * (i.e. This includes ONLY EXONS.)
     */
    private String  alignedCodingSequenceReferenceAllele = null;

    /**
     * The in-frame position (1-based, inclusive) of the last base of the last codon that contains the reference
     * allele relative to the start of the coding sequence.
     * All codons containing the reference allele can be extracted from the coding sequence using
     * {@link SequenceComparison#alignedCodingSequenceAlleleStart} and {@link SequenceComparison#alignedReferenceAlleleStop}.
     */
    private Integer alignedReferenceAlleleStop           = null;

    /**
     * The {@link String} representation of the alternate allele, correctly represented for the {@link Strand} on which it appears.
     */
    private String  alternateAllele                      = null;

    /**
     * An in-frame sequence of bases that includes the entire alternate allele based on the reference genome.
     * May span multiple codons.
     * May include intron bases.
     */
    private String  alignedAlternateAllele               = null;

    /**
     * An in-frame sequence of bases that includes the entire alternate allele based on the coding sequence reference.
     * May span multiple codons.
     */
    private String  alignedCodingSequenceAlternateAllele = null;

    /**
     * The in-frame position (1-based, inclusive) of the last base of the last codon that contains the alternate
     * allele relative to the start of the coding sequence.
     * All codons containing the alternate allele can be extracted from the coding sequence using
     * {@link SequenceComparison#alignedCodingSequenceAlleleStart} and {@link SequenceComparison#alignedAlternateAlleleStop}.
     */
    private Integer alignedAlternateAlleleStop           = null;

    /**
     * The fraction of Guanine and Cytosine bases in a window of a given size around a variant.
     * The default windows size is {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory#gcContentWindowSizeBases}.
     */
    private Double gcContent                             = null;

    /**
     * Container holding information on the protein change represented by this {@link SequenceComparison} object.
     */
    private ProteinChangeInfo proteinChangeInfo = null;

    // =============================================================================================================

    public boolean hasSequenceInfo() {
        return this.transcriptCodingSequence != null;
    }

    public String getReferenceBases() {
        return referenceBases;
    }

    public SequenceComparison setReferenceBases(final String referenceBases) {
        this.referenceBases = referenceBases;
        return this;
    }

    public SequenceComparison setReferenceWindow(final Integer referenceWindow) {
        this.referenceWindow = referenceWindow;
        return this;
    }

   /**
     * Return the {@link ReferenceSequence} containing the coding region for the transcript of this {@link SequenceComparison}.
     * This does NOT include introns.
     * The reference sequence is stored in the forward reading direction.
     * For NEGATIVE strand reads, must reverse complement any bases retrieved.
     */
    public ReferenceSequence getTranscriptCodingSequence() {
        return transcriptCodingSequence;
    }

     public SequenceComparison setTranscriptCodingSequence(final ReferenceSequence transcriptCodingSequence) {
        this.transcriptCodingSequence = transcriptCodingSequence;
         return this;
    }

    public String getContig() {
        return contig;
    }

    public SequenceComparison setContig(final String contig) {
        this.contig = contig;
        return this;
    }

    public Strand getStrand() {
        return strand;
    }

    public SequenceComparison setStrand(final Strand strand) {

        if (strand == Strand.NONE) {
            throw new GATKException("Cannot handle NONE strand.");
        }

        this.strand = strand;

        return this;
    }

    public Integer getAlleleStart() {
        return alleleStart;
    }

    public SequenceComparison setAlleleStart(final Integer alleleStart) {
        this.alleleStart = alleleStart;
        return this;
    }

    public Integer getTranscriptAlleleStart() {
        return transcriptAlleleStart;
    }

    public SequenceComparison setTranscriptAlleleStart(final Integer transcriptAlleleStart) {
        this.transcriptAlleleStart = transcriptAlleleStart;
        return this;
    }

    public Integer getCodingSequenceAlleleStart() {
        return codingSequenceAlleleStart;
    }

    public SequenceComparison setCodingSequenceAlleleStart(final Integer codingSequenceAlleleStart) {
        this.codingSequenceAlleleStart = codingSequenceAlleleStart;
        return this;
    }

    public Integer getAlignedCodingSequenceAlleleStart() {
        return alignedCodingSequenceAlleleStart;
    }

    public SequenceComparison setAlignedCodingSequenceAlleleStart(final Integer alignedCodingSequenceAlleleStart) {
        this.alignedCodingSequenceAlleleStart = alignedCodingSequenceAlleleStart;
        return this;
    }

    public Integer getExonStartPosition() {
        return exonStartPosition;
    }

    public Integer getExonEndPosition() {
        return exonEndPosition;
    }

    public SequenceComparison setExonPosition( final SimpleInterval exonPosition ) {
        this.exonStartPosition = exonPosition.getStart();
        this.exonEndPosition = exonPosition.getEnd();

        return this;
    }

    public String getReferenceAllele() {
        return referenceAllele;
    }

    public SequenceComparison setReferenceAllele(final String referenceAllele) {
        this.referenceAllele = referenceAllele;
        return this;
    }

    public String getAlignedReferenceAllele() {
        return alignedReferenceAllele;
    }

    public SequenceComparison setAlignedReferenceAllele(final String alignedReferenceAllele) {
        this.alignedReferenceAllele = alignedReferenceAllele;
        return this;
    }

    public String getAlignedCodingSequenceReferenceAllele() {
        return alignedCodingSequenceReferenceAllele;
    }

    public SequenceComparison setAlignedCodingSequenceReferenceAllele(final String alignedCodingSequenceReferenceAllele) {
        this.alignedCodingSequenceReferenceAllele = alignedCodingSequenceReferenceAllele;
        return this;
    }

    public Integer getAlignedReferenceAlleleStop() {
        return alignedReferenceAlleleStop;
    }

    public SequenceComparison setAlignedReferenceAlleleStop(final Integer alignedReferenceAlleleStop) {
        this.alignedReferenceAlleleStop = alignedReferenceAlleleStop;
        return this;
    }

    public String getAlternateAllele() {
        return alternateAllele;
    }

    public SequenceComparison setAlternateAllele(final String alternateAllele) {
        this.alternateAllele = alternateAllele;
        return this;
    }

    public String getAlignedAlternateAllele() {
        return alignedAlternateAllele;
    }

    public SequenceComparison setAlignedAlternateAllele(final String alignedAlternateAllele) {
        this.alignedAlternateAllele = alignedAlternateAllele;
        return this;
    }

    public String getAlignedCodingSequenceAlternateAllele() {
        return alignedCodingSequenceAlternateAllele;
    }

    public SequenceComparison setAlignedCodingSequenceAlternateAllele(final String alignedCodingSequenceAlternateAllele) {
        this.alignedCodingSequenceAlternateAllele = alignedCodingSequenceAlternateAllele;
        return this;
    }

    public SequenceComparison setAlignedAlternateAlleleStop(final Integer alignedAlternateAlleleStop) {
        this.alignedAlternateAlleleStop = alignedAlternateAlleleStop;
        return this;
    }

    public Double getGcContent() {
        return gcContent;
    }

    public SequenceComparison setGcContent(final Double gcContent) {
        this.gcContent = gcContent;
        return this;
    }

    public SequenceComparison setProteinChangeInfo(final ProteinChangeInfo p) { proteinChangeInfo = p; return this;}

    public ProteinChangeInfo getProteinChangeInfo() {
        return proteinChangeInfo;
    }
}