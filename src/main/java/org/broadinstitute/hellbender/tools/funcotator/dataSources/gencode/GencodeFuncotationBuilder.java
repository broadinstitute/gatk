package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfGeneFeature;

import java.util.ArrayList;
import java.util.List;

/**
 * A builder object to create {@link GencodeFuncotation}s.
 * This follows the Builder Pattern, and as such is a low-level object that is separate from {@link GencodeFuncotationFactory}.
 * Created by jonn on 11/6/17.
 */
public class GencodeFuncotationBuilder {

    //==================================================================================================================

    /**
     * The {@link GencodeFuncotation} that is being populated by this {@link GencodeFuncotationBuilder}.
     */
    private GencodeFuncotation gencodeFuncotation;

    //==================================================================================================================

    public GencodeFuncotationBuilder() {
        gencodeFuncotation = new GencodeFuncotation();
    }

    public GencodeFuncotationBuilder(final GencodeFuncotation gf) {
        gencodeFuncotation = new GencodeFuncotation(gf);
    }

    public GencodeFuncotation build() {
        // TODO: In the future, we will need a mechanism for populating the metadata  (https://github.com/broadinstitute/gatk/issues/4857)
        gencodeFuncotation.setMetadata(FuncotationMetadataUtils.createWithUnknownAttributes(new ArrayList<>(gencodeFuncotation.getFieldNames())));
        return gencodeFuncotation;
    }

    //==================================================================================================================

    /**
     * Set the Reference Allele in the {@link GencodeFuncotation}.
     * @param refAllele The reference {@link Allele} for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setRefAllele( final Allele refAllele ) {
        gencodeFuncotation.setRefAllele(refAllele.getBaseString());
        return this;
    }

    /**
     * Set the Strand in the {@link GencodeFuncotation}.
     * @param strand The reference {@link Strand} for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setStrand( final Strand strand ) {
        switch(strand) {
            case POSITIVE: gencodeFuncotation.setTranscriptStrand( "+" ); break;
            case NEGATIVE: gencodeFuncotation.setTranscriptStrand( "-"); break;
            default:
                throw new GATKException.ShouldNeverReachHereException("An invalid Strand type was specified.  Strand must be POSITIVE or NEGATIVE.");
        }
        return this;
    }

    /**
     * Set the Hugo Symbol / Gene Name in the {@link GencodeFuncotation}.
     * @param hugoSymbol The Hugo Symbol / Gene Name for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setHugoSymbol(final String hugoSymbol) {
        gencodeFuncotation.setHugoSymbol(hugoSymbol);
        return this;
    }

    /**
     * Set the NCBI build in the {@link GencodeFuncotation}.
     * @param ncbiBuild The NCBI build for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setNcbiBuild(final String ncbiBuild) {
        gencodeFuncotation.setNcbiBuild(ncbiBuild);
        return this;
    }

    /**
     * Set the Chromosome name in the {@link GencodeFuncotation}.
     * @param chromosomeName The Chromosome name for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setChromosome(final String chromosomeName) {
        gencodeFuncotation.setChromosome(chromosomeName);
        return this;
    }

    /**
     * Set the start position (1-based, inclusive) in the {@link GencodeFuncotation}.
     * @param start The start position (1-based, inclusive) for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setStart(final int start) {
        gencodeFuncotation.setStart(start);
        return this;
    }

    /**
     * Set the end position (1-based, inclusive) in the {@link GencodeFuncotation}.
     * @param end The end position (1-based, inclusive) for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setEnd(final int end) {
        gencodeFuncotation.setEnd(end);
        return this;
    }

    /**
     * Set the Variant Type (1-based, inclusive) in the {@link GencodeFuncotation}.
     * @param variantType The {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantType} for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setVariantType(final GencodeFuncotation.VariantType variantType) {
        gencodeFuncotation.setVariantType(variantType);
        return this;
    }

    /**
     * Set the bases of tumorSeqAllele2 in the {@link GencodeFuncotation}.
     * @param altAllele A {@link String} containing the bases in the allele to set as tumorSeqAllele2 for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setTumorSeqAllele2(final String altAllele) {
        gencodeFuncotation.setTumorSeqAllele2(altAllele);
        return this;
    }

    /**
     * Set the Genome Change in the {@link GencodeFuncotation}.
     * @param genomeChange The Genome Change for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setGenomeChange(final String genomeChange) {
        gencodeFuncotation.setGenomeChange(genomeChange);
        return this;
    }

    /**
     * Set the Transcript ID in the {@link GencodeFuncotation}.
     * @param transcriptId The {@link String} containing the transcript ID for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setAnnotationTranscript(final String transcriptId) {
        gencodeFuncotation.setAnnotationTranscript(transcriptId);
        return this;
    }

    /**
     * Set the position (1-based, inclusive) relative to the start of the transcript of a the variant in the {@link GencodeFuncotation}.
     * @param transcriptPos The position (1-based, inclusive) relative to the start of the transcript of a the variant in the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setTranscriptPos(final Integer transcriptPos) {
        gencodeFuncotation.setTranscriptPos(transcriptPos);
        return this;
    }

    /**
     * Set the list of other Transcript IDs in the {@link GencodeFuncotation}.
     * @param transcripts The {@link List} of {@link String}s containing the transcript IDs for the other transcripts for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setOtherTranscripts(final List<String> transcripts) {
        gencodeFuncotation.setOtherTranscripts(transcripts);
        return this;
    }

    /**
     * Set the Variant Classification in the {@link GencodeFuncotation}.
     * @param variantClassification The {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setVariantClassification(final GencodeFuncotation.VariantClassification variantClassification) {
        gencodeFuncotation.setVariantClassification( variantClassification );
        return this;
    }

    /**
     * Set the Secondary Variant Classification in the {@link GencodeFuncotation}.
     * @param variantClassification The Secondary {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setSecondaryVariantClassification(final GencodeFuncotation.VariantClassification variantClassification) {
        gencodeFuncotation.setSecondaryVariantClassification( variantClassification );
        return this;
    }

    /**
     * Set the transcript exon number in the {@link GencodeFuncotation}.
     * @param transcriptExonNumber The number (>=1) of the exon in the transcript for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setTranscriptExonNumber( final int transcriptExonNumber ) {
        gencodeFuncotation.setTranscriptExonNumber( transcriptExonNumber );
        return this;
    }

    /**
     * Set the Codon Change {@link String} in the {@link GencodeFuncotation}.
     * @param codonChange The {@link String} containing the Codon Change information for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setCodonChange( final String codonChange ) {
        gencodeFuncotation.setCodonChange( codonChange );
        return this;
    }

    /**
     * Set the Protein Change {@link String} in the {@link GencodeFuncotation}.
     * @param proteinChange The {@link String} containing the Protein Change information for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setProteinChange( final String proteinChange ) {
        gencodeFuncotation.setProteinChange( proteinChange );
        return this;
    }

    /**
     * Set the CDNA Change {@link String} in the {@link GencodeFuncotation}.
     * @param cDnaChange The {@link String} containing the CDNA Change information for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setcDnaChange( final String cDnaChange ) {
        gencodeFuncotation.setcDnaChange( cDnaChange );
        return this;
    }

    /**
     * Set the LocusLevel {@link Integer} in the {@link GencodeFuncotation}.
     * @param locusLevel The {@link Integer} containing the locus level for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setLocusLevel( final Integer locusLevel ) {
        gencodeFuncotation.setLocusLevel( locusLevel );
        return this;
    }

    /**
     * Set the Appris Rank in the {@link GencodeFuncotation}.
     * @param apprisRank The {@link GencodeGtfGeneFeature.FeatureTag} containing the Appris Rank for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setApprisRank( final GencodeGtfGeneFeature.FeatureTag apprisRank ) {
        gencodeFuncotation.setApprisRank( apprisRank );
        return this;
    }

    /**
     * Set the Transcript Length {@link Integer} in the {@link GencodeFuncotation}.
     * @param transcriptLength The {@link Integer} containing the Transcript Length for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setTranscriptLength( final Integer transcriptLength ) {
        gencodeFuncotation.setTranscriptLength( transcriptLength );
        return this;
    }

    /**
     * Set the GC Content {@link Double} in the {@link GencodeFuncotation}.
     * @param gcContent The {@link Double} containing the GC Content for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setGcContent( final Double gcContent ) {
        gencodeFuncotation.setGcContent( gcContent );
        return this;
    }

    /**
     * Set the Reference Context {@link String} in the {@link GencodeFuncotation}.
     * @param referenceContext The {@link String} containing the ReferenceContext for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setReferenceContext(final String referenceContext) {
        gencodeFuncotation.setReferenceContext( referenceContext );
        return this;
    }

    /**
     * Set the Version {@link String} in the {@link GencodeFuncotation}.
     * @param version The {@link String} containing the Version for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setVersion(final String version) {
        gencodeFuncotation.setVersion( version );
        return this;
    }

    /**
     * Set the GeneTranscriptType {@link org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature.GeneTranscriptType} in the {@link GencodeFuncotation}.
     * @param geneTranscriptType The {@link String} containing the GeneTranscriptType for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setGeneTranscriptType(final GencodeGtfFeature.GeneTranscriptType geneTranscriptType) {
        gencodeFuncotation.setGeneTranscriptType( geneTranscriptType );
        return this;
    }

    /**
     * Set the Data Source Name {@link String} in the {@link GencodeFuncotation}.
     * @param name The {@link String} containing the Data Source Name for the {@link GencodeFuncotation}.
     * @return {@code this} {@link GencodeFuncotationBuilder}
     */
    public GencodeFuncotationBuilder setDataSourceName( final String name ) {
        gencodeFuncotation.setDataSourceName( name );
        return this;
    }
 }
