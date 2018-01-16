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
     * The position (1-based, inclusive in genome coordinates - relative to the start of
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
     * The start position (1-based, inclusive) of the Protein Change for the alleles in this {@link SequenceComparison}.
     * This is computed using {@link SequenceComparison#alignedCodingSequenceAlleleStart}.
     * See {@link FuncotatorUtils#getProteinChangePosition} for more details.
     */
    private Integer proteinChangeStartPosition           = null;

    /**
     * The end position (1-based, inclusive) of the Protein Change for the alleles in this {@link SequenceComparison}.
     * This is computed by using {@link SequenceComparison#alignedCodingSequenceAlleleStart} and the length of
     * {@link SequenceComparison#alignedCodingSequenceAlternateAllele}.
     * See {@link FuncotatorUtils#getProteinChangeEndPosition} for more details.
     */
    private Integer proteinChangeEndPosition             = null;

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
     * The amino acid sequence as coded by the in-frame reference coding sequence ({@link SequenceComparison#alignedCodingSequenceReferenceAllele}.
     */
    private String  referenceAminoAcidSequence           = null;

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
     * The amino acid sequence as coded by the in-frame alternate coding sequence ({@link SequenceComparison#alignedCodingSequenceAlternateAllele}.
     */
    private String  alternateAminoAcidSequence           = null;

    /**
     * The fraction of Guanine and Cytosine bases in a window of a given size around a variant.
     * The default windows size is {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory#gcContentWindowSizeBases}.
     */
    private Double gcContent                             = null;

    // =============================================================================================================

    public boolean hasSequenceInfo() {
        return this.transcriptCodingSequence != null;
    }

    public String getReferenceBases() {
        return referenceBases;
    }

    public void setReferenceBases(final String referenceBases) {
        this.referenceBases = referenceBases;
    }

    public Integer getReferenceWindow() {
        return referenceWindow;
    }

    public void setReferenceWindow(final Integer referenceWindow) {
        this.referenceWindow = referenceWindow;
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

     public void setTranscriptCodingSequence(final ReferenceSequence transcriptCodingSequence) {
        this.transcriptCodingSequence = transcriptCodingSequence;
    }

    /**
     * Return the {@link ReferenceSequence} containing the coding region for the transcript of this {@link SequenceComparison}.
     * This does NOT include introns.
     * The reference sequence is stored in the forward reading direction.
     * For NEGATIVE strand reads, must reverse complement any bases retrieved.
     */
    public ReferenceSequence getReferenceCodingSequence() {
        return referenceCodingSequence;
    }

    public void setReferenceCodingSequence(final ReferenceSequence referenceCodingSequence) {
        this.referenceCodingSequence = referenceCodingSequence;
    }

    public String getContig() {
        return contig;
    }

    public void setContig(final String contig) {
        this.contig = contig;
    }

    public Strand getStrand() {
        return strand;
    }

    public void setStrand(final Strand strand) {

        if (strand == Strand.NONE) {
            throw new GATKException("Cannot handle NONE strand.");
        }

        this.strand = strand;
    }

    public Integer getAlleleStart() {
        return alleleStart;
    }

    public void setAlleleStart(final Integer alleleStart) {
        this.alleleStart = alleleStart;
    }

    public Integer getTranscriptAlleleStart() {
        return transcriptAlleleStart;
    }

    public void setTranscriptAlleleStart(final Integer transcriptAlleleStart) {
        this.transcriptAlleleStart = transcriptAlleleStart;
    }

    public Integer getCodingSequenceAlleleStart() {
        return codingSequenceAlleleStart;
    }

    public void setCodingSequenceAlleleStart(final Integer codingSequenceAlleleStart) {
        this.codingSequenceAlleleStart = codingSequenceAlleleStart;
    }

    public Integer getAlignedCodingSequenceAlleleStart() {
        return alignedCodingSequenceAlleleStart;
    }

    public void setAlignedCodingSequenceAlleleStart(final Integer alignedCodingSequenceAlleleStart) {
        this.alignedCodingSequenceAlleleStart = alignedCodingSequenceAlleleStart;
    }

    public Integer getExonStartPosition() {
        return exonStartPosition;
    }

    public void setExonStartPosition(final Integer exonStartPosition) {
        this.exonStartPosition = exonStartPosition;
    }

    public Integer getExonEndPosition() {
        return exonEndPosition;
    }

    public void setExonEndPosition(final Integer exonEndPosition) {
        this.exonEndPosition = exonEndPosition;
    }

    public void setExonPosition( final SimpleInterval exonPosition ) {
        this.exonStartPosition = exonPosition.getStart();
        this.exonEndPosition = exonPosition.getEnd();
    }

    public Integer getProteinChangeStartPosition() {
        return proteinChangeStartPosition;
    }

    public void setProteinChangeStartPosition(final Integer proteinChangeStartPosition) {
        this.proteinChangeStartPosition = proteinChangeStartPosition;
    }

    public Integer getProteinChangeEndPosition() {
        return proteinChangeEndPosition;
    }

    public void setProteinChangeEndPosition(final Integer proteinChangeEndPosition) {
        this.proteinChangeEndPosition = proteinChangeEndPosition;
    }

    public String getReferenceAllele() {
        return referenceAllele;
    }

    public void setReferenceAllele(final String referenceAllele) {
        this.referenceAllele = referenceAllele;
    }

    public String getAlignedReferenceAllele() {
        return alignedReferenceAllele;
    }

    public void setAlignedReferenceAllele(final String alignedReferenceAllele) {
        this.alignedReferenceAllele = alignedReferenceAllele;
    }

    public String getAlignedCodingSequenceReferenceAllele() {
        return alignedCodingSequenceReferenceAllele;
    }

    public void setAlignedCodingSequenceReferenceAllele(final String alignedCodingSequenceReferenceAllele) {
        this.alignedCodingSequenceReferenceAllele = alignedCodingSequenceReferenceAllele;
    }

    public Integer getAlignedReferenceAlleleStop() {
        return alignedReferenceAlleleStop;
    }

    public void setAlignedReferenceAlleleStop(final Integer alignedReferenceAlleleStop) {
        this.alignedReferenceAlleleStop = alignedReferenceAlleleStop;
    }

    public String getReferenceAminoAcidSequence() {
        return referenceAminoAcidSequence;
    }

    public void setReferenceAminoAcidSequence(final String referenceAminoAcidSequence) {
        this.referenceAminoAcidSequence = referenceAminoAcidSequence;
    }

    public String getAlternateAllele() {
        return alternateAllele;
    }

    public void setAlternateAllele(final String alternateAllele) {
        this.alternateAllele = alternateAllele;
    }

    public String getAlignedAlternateAllele() {
        return alignedAlternateAllele;
    }

    public void setAlignedAlternateAllele(final String alignedAlternateAllele) {
        this.alignedAlternateAllele = alignedAlternateAllele;
    }

    public String getAlignedCodingSequenceAlternateAllele() {
        return alignedCodingSequenceAlternateAllele;
    }

    public void setAlignedCodingSequenceAlternateAllele(final String alignedCodingSequenceAlternateAllele) {
        this.alignedCodingSequenceAlternateAllele = alignedCodingSequenceAlternateAllele;
    }

    public Integer getAlignedAlternateAlleleStop() {
        return alignedAlternateAlleleStop;
    }

    public void setAlignedAlternateAlleleStop(final Integer alignedAlternateAlleleStop) {
        this.alignedAlternateAlleleStop = alignedAlternateAlleleStop;
    }

    public String getAlternateAminoAcidSequence() {
        return alternateAminoAcidSequence;
    }

    public void setAlternateAminoAcidSequence(final String alternateAminoAcidSequence) {
        this.alternateAminoAcidSequence = alternateAminoAcidSequence;
    }

    public Double getGcContent() {
        return gcContent;
    }

    public void setGcContent(final Double gcContent) {
        this.gcContent = gcContent;
    }
}