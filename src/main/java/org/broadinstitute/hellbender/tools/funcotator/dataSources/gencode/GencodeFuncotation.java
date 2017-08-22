package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import org.broadinstitute.hellbender.tools.funcotator.Funcotation;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A class to represent a Functional Annotation.
 * Created by jonn on 8/22/17.
 */
public class GencodeFuncotation extends Funcotation {

    private static final String FIELD_DELIMITER = "|";
    private static final String OTHER_TRANSCRIPT_DELIMITER = ";";

    //==================================================================================================================

    private String                  hugoSymbol;                         // TRIVIAL (i.e. by the time we match to a transcript, we have this info regardless to where in the transcript the variant lies)
    private String                  ncbiBuild;                          // TRIVIAL
    private String                  chromosome;                         // TRIVIAL
    private int                     start;                              // TRIVIAL
    private int                     end;                                // TRIVIAL
    private VariantClassification   variantClassification;              //          CDS / UTR / INTRON / IGR
    private VariantClassification   secondaryVariantClassification;     //          CDS / INTRON
    private VariantType             variantType;                        // TRIVIAL
    private String                  refAllele;                          // TRIVIAL
    private String                  tumorSeqAllele1;                    // TRIVIAL
    private String                  tumorSeqAllele2;                    // TRIVIAL

    private String                  genomeChange;                       // TRIVIAL
    private String                  annotationTranscript;               // TRIVIAL
    private String                  transcriptStrand;                   // TRIVIAL
    private Integer                 transcriptExon;                     //           CDS / UTRs
    private Integer                 transcriptPos;                      // TRIVIAL
    private String                  cDnaChange;                         //           CDS
    private String                  codonChange;                        //           CDS
    private String                  proteinChange;                      //           CDS
    private List<String>            otherTranscripts;                   // TRIVIAL

    //==================================================================================================================

    /**
     * Basic constructor for a {@link GencodeFuncotation}.
     */
    public GencodeFuncotation() {}

    //==================================================================================================================

    /**
     * @return An ordered {@link List} of {@link String} containing the field names that {@link GencodeFuncotation} produces.
     */
    public static List<String> getSerializedFieldNames() {

        final List<String> fields = new ArrayList<>();

        for(final Field f : GencodeFuncotation.class.getDeclaredFields() ) {
            if ( !Modifier.isStatic(f.getModifiers()) ) {
                fields.add( f.getName() );
            }
        }

        return fields;
    }

    //==================================================================================================================

    /**
     * Converts this {@link GencodeFuncotation} to a string suitable for insertion into a VCF file.
     * @return a {@link String} representing this {@link GencodeFuncotation} suitable for insertion into a VCF file.
     */
    @Override
    public String serializeToVcfString() {

        return (hugoSymbol != null ? hugoSymbol : "") + FIELD_DELIMITER +
                (ncbiBuild != null ? ncbiBuild : "") + FIELD_DELIMITER +
                (chromosome != null ? chromosome : "") + FIELD_DELIMITER +
                start + FIELD_DELIMITER +
                end + FIELD_DELIMITER +
                (variantClassification != null ? variantClassification : "") + FIELD_DELIMITER +
                (secondaryVariantClassification != null ? secondaryVariantClassification : "") + FIELD_DELIMITER +
                (variantType != null ? variantType : "") + FIELD_DELIMITER +
                (refAllele != null ? refAllele : "") + FIELD_DELIMITER +
                (tumorSeqAllele1 != null ? tumorSeqAllele1 : "") + FIELD_DELIMITER +
                (tumorSeqAllele2 != null ? tumorSeqAllele2 : "") + FIELD_DELIMITER +
                (genomeChange != null ? genomeChange : "") + FIELD_DELIMITER +
                (annotationTranscript != null ? annotationTranscript : "") + FIELD_DELIMITER +
                (transcriptStrand != null ? transcriptStrand : "") + FIELD_DELIMITER +
                (transcriptExon != null ? transcriptExon : "") + FIELD_DELIMITER +
                (transcriptPos != null ? transcriptPos : "") + FIELD_DELIMITER +
                (cDnaChange != null ? cDnaChange : "") + FIELD_DELIMITER +
                (codonChange != null ? codonChange : "") + FIELD_DELIMITER +
                (proteinChange != null ? proteinChange : "") + FIELD_DELIMITER +
                (otherTranscripts != null ? otherTranscripts.stream().map(Object::toString).collect(Collectors.joining(OTHER_TRANSCRIPT_DELIMITER)) : "");
    }

    //==================================================================================================================

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final GencodeFuncotation that = (GencodeFuncotation) o;

        if (start != that.start) return false;
        if (end != that.end) return false;
        if (hugoSymbol != null ? !hugoSymbol.equals(that.hugoSymbol) : that.hugoSymbol != null) return false;
        if (ncbiBuild != null ? !ncbiBuild.equals(that.ncbiBuild) : that.ncbiBuild != null) return false;
        if (chromosome != null ? !chromosome.equals(that.chromosome) : that.chromosome != null) return false;
        if (variantClassification != that.variantClassification) return false;
        if (secondaryVariantClassification != that.secondaryVariantClassification) return false;
        if (variantType != that.variantType) return false;
        if (refAllele != null ? !refAllele.equals(that.refAllele) : that.refAllele != null) return false;
        if (tumorSeqAllele1 != null ? !tumorSeqAllele1.equals(that.tumorSeqAllele1) : that.tumorSeqAllele1 != null)
            return false;
        if (tumorSeqAllele2 != null ? !tumorSeqAllele2.equals(that.tumorSeqAllele2) : that.tumorSeqAllele2 != null)
            return false;
        if (genomeChange != null ? !genomeChange.equals(that.genomeChange) : that.genomeChange != null) return false;
        if (annotationTranscript != null ? !annotationTranscript.equals(that.annotationTranscript) : that.annotationTranscript != null)
            return false;
        if (transcriptStrand != null ? !transcriptStrand.equals(that.transcriptStrand) : that.transcriptStrand != null)
            return false;
        if (transcriptExon != null ? !transcriptExon.equals(that.transcriptExon) : that.transcriptExon != null)
            return false;
        if (transcriptPos != null ? !transcriptPos.equals(that.transcriptPos) : that.transcriptPos != null)
            return false;
        if (cDnaChange != null ? !cDnaChange.equals(that.cDnaChange) : that.cDnaChange != null) return false;
        if (codonChange != null ? !codonChange.equals(that.codonChange) : that.codonChange != null) return false;
        if (proteinChange != null ? !proteinChange.equals(that.proteinChange) : that.proteinChange != null)
            return false;
        return otherTranscripts != null ? otherTranscripts.equals(that.otherTranscripts) : that.otherTranscripts == null;
    }

    @Override
    public int hashCode() {
        int result = hugoSymbol != null ? hugoSymbol.hashCode() : 0;
        result = 31 * result + (ncbiBuild != null ? ncbiBuild.hashCode() : 0);
        result = 31 * result + (chromosome != null ? chromosome.hashCode() : 0);
        result = 31 * result + start;
        result = 31 * result + end;
        result = 31 * result + (variantClassification != null ? variantClassification.hashCode() : 0);
        result = 31 * result + (secondaryVariantClassification != null ? secondaryVariantClassification.hashCode() : 0);
        result = 31 * result + (variantType != null ? variantType.hashCode() : 0);
        result = 31 * result + (refAllele != null ? refAllele.hashCode() : 0);
        result = 31 * result + (tumorSeqAllele1 != null ? tumorSeqAllele1.hashCode() : 0);
        result = 31 * result + (tumorSeqAllele2 != null ? tumorSeqAllele2.hashCode() : 0);
        result = 31 * result + (genomeChange != null ? genomeChange.hashCode() : 0);
        result = 31 * result + (annotationTranscript != null ? annotationTranscript.hashCode() : 0);
        result = 31 * result + (transcriptStrand != null ? transcriptStrand.hashCode() : 0);
        result = 31 * result + (transcriptExon != null ? transcriptExon.hashCode() : 0);
        result = 31 * result + (transcriptPos != null ? transcriptPos.hashCode() : 0);
        result = 31 * result + (cDnaChange != null ? cDnaChange.hashCode() : 0);
        result = 31 * result + (codonChange != null ? codonChange.hashCode() : 0);
        result = 31 * result + (proteinChange != null ? proteinChange.hashCode() : 0);
        result = 31 * result + (otherTranscripts != null ? otherTranscripts.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "GencodeFuncotation{" +
                "hugoSymbol='" + hugoSymbol + '\'' +
                ", ncbiBuild='" + ncbiBuild + '\'' +
                ", chromosome='" + chromosome + '\'' +
                ", start=" + start +
                ", end=" + end +
                ", variantClassification=" + variantClassification +
                ", secondaryVariantClassification=" + secondaryVariantClassification +
                ", variantType=" + variantType +
                ", refAllele='" + refAllele + '\'' +
                ", tumorSeqAllele1='" + tumorSeqAllele1 + '\'' +
                ", tumorSeqAllele2='" + tumorSeqAllele2 + '\'' +
                ", genomeChange='" + genomeChange + '\'' +
                ", annotationTranscript='" + annotationTranscript + '\'' +
                ", transcriptStrand='" + transcriptStrand + '\'' +
                ", transcriptExon=" + transcriptExon +
                ", transcriptPos=" + transcriptPos +
                ", cDnaChange='" + cDnaChange + '\'' +
                ", codonChange='" + codonChange + '\'' +
                ", proteinChange='" + proteinChange + '\'' +
                ", otherTranscripts=" + otherTranscripts +
                '}';
    }

//==================================================================================================================

    public String getHugoSymbol() {
        return hugoSymbol;
    }

    public void setHugoSymbol(final String hugoSymbol) {
        this.hugoSymbol = hugoSymbol;
    }

    public String getNcbiBuild() {
        return ncbiBuild;
    }

    public void setNcbiBuild(final String ncbiBuild) {
        this.ncbiBuild = ncbiBuild;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(final String chromosome) {
        this.chromosome = chromosome;
    }

    public int getStart() {
        return start;
    }

    public void setStart(final int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(final int end) {
        this.end = end;
    }

    public VariantClassification getVariantClassification() {
        return variantClassification;
    }

    public void setVariantClassification(final VariantClassification variantClassification) {
        this.variantClassification = variantClassification;
    }

    public VariantClassification getSecondaryVariantClassification() {
        return secondaryVariantClassification;
    }

    public void setSecondaryVariantClassification(final VariantClassification secondaryVariantClassification) {
        this.secondaryVariantClassification = secondaryVariantClassification;
    }

    public VariantType getVariantType() {
        return variantType;
    }

    public void setVariantType(final VariantType variantType) {
        this.variantType = variantType;
    }

    public String getRefAllele() {
        return refAllele;
    }

    public void setRefAllele(final String refAllele) {
        this.refAllele = refAllele;
    }

    public String getTumorSeqAllele1() {
        return tumorSeqAllele1;
    }

    public void setTumorSeqAllele1(final String tumorSeqAllele1) {
        this.tumorSeqAllele1 = tumorSeqAllele1;
    }

    public String getTumorSeqAllele2() {
        return tumorSeqAllele2;
    }

    public void setTumorSeqAllele2(final String tumorSeqAllele2) {
        this.tumorSeqAllele2 = tumorSeqAllele2;
    }

    public String getGenomeChange() {
        return genomeChange;
    }

    public void setGenomeChange(final String genomeChange) {
        this.genomeChange = genomeChange;
    }

    public String getAnnotationTranscript() {
        return annotationTranscript;
    }

    public void setAnnotationTranscript(final String annotationTranscript) {
        this.annotationTranscript = annotationTranscript;
    }

    public String getTranscriptStrand() {
        return transcriptStrand;
    }

    public void setTranscriptStrand(final String transcriptStrand) {
        this.transcriptStrand = transcriptStrand;
    }

    public Integer getTranscriptExon() {
        return transcriptExon;
    }

    public void setTranscriptExon(final Integer transcriptExon) {
        this.transcriptExon = transcriptExon;
    }

    public Integer getTranscriptPos() {
        return transcriptPos;
    }

    public void setTranscriptPos(final Integer transcriptPos) {
        this.transcriptPos = transcriptPos;
    }

    public String getcDnaChange() {
        return cDnaChange;
    }

    public void setcDnaChange(final String cDnaChange) {
        this.cDnaChange = cDnaChange;
    }

    public String getCodonChange() {
        return codonChange;
    }

    public void setCodonChange(final String codonChange) {
        this.codonChange = codonChange;
    }

    public String getProteinChange() {
        return proteinChange;
    }

    public void setProteinChange(final String proteinChange) {
        this.proteinChange = proteinChange;
    }

    public List<String> getOtherTranscripts() {
        return otherTranscripts;
    }

    public void setOtherTranscripts(final List<String> otherTranscripts) {
        this.otherTranscripts = otherTranscripts;
    }


    //==================================================================================================================

    public enum VariantType {
        INS("INS"),
        DEL("DEL"),
        SNP("SNP"),
        DNP("DNP"),
        TNP("TNP"),
        ONP("ONP"),
        MNP("MNP"),
        xNP("NP");

        final private String serialized;

        VariantType(final String serialized) { this.serialized = serialized; }

        @Override
        public String toString() {
            return serialized;
        }
    }

    /**
     * Represents the type and severity of a variant.
     * The lower the {@link VariantClassification#relativeSeverity}, the more severe the variant.
     * Descriptions taken from:
     *     https://gatkforums.broadinstitute.org/gatk/discussion/8815/oncotator-variant-classification-and-secondary-variant-classification
     */
    public enum VariantClassification {
        // Only for Introns:
        /** Variant lies between exons within the bounds of the chosen transcript. */
        INTRON(10),

        // Only for UTRs:
        /** Variant is on the 5'UTR for the chosen transcript. */
        FIVE_PRIME_UTR(6),
        /** Variant is on the 3'UTR for the chosen transcript */
        THREE_PRIME_UTR(6),

        // Only for IGRs:
        /** Intergenic region. Does not overlap any transcript. */
        IGR(20),
        /** The variant is upstream of the chosen transcript (within 3kb) */
        FIVE_PRIME_FLANK(15),

        // These can be in Coding regions or Introns (only SPLICE_SITE):
        /** The point mutation alters the protein structure by one amino acid. */
        MISSENSE(1),
        /** A premature stop codon is created by the variant. */
        NONSENSE(0),
        /** Variant removes stop codon. */
        NONSTOP(0),
        /** Variant is in coding region of the chosen transcript, but protein structure is identical (i.e. a synonymous mutation). */
        SILENT(5),
        /** The variant is within a configurable number of bases (default=2) of a splice site. See the secondary classification to determine if it lies on the exon or intron side. */
        SPLICE_SITE(4),
        /** Deletion that keeps the sequence in frame (i.e. deletion of a length evenly divisible by 3). */
        IN_FRAME_DEL(1),
        /** Insertion that keeps the sequence in frame (i.e. insertion of a length evenly divisible by 3). */
        IN_FRAME_INS(1),
        /** Insertion that moves the coding sequence out of frame (i.e. insertion of a length not evenly divisible by 3). */
        FRAME_SHIFT_INS(2),
        /** Deletion that moves the sequence out of frame (i.e. deletion of a length not evenly divisible by 3). */
        FRAME_SHIFT_DEL(2),
        /** Point mutation that overlaps the start codon. */
        START_CODON_SNP(3),
        /** Insertion that overlaps the start codon. */
        START_CODON_INS(3),
        /** Deletion that overlaps the start codon. */
        START_CODON_DEL(3),

        // These can only be in 5' UTRs:
        /** New start codon is created by the given variant using the chosen transcript. However, it is in frame relative to the coded protein. */
        DE_NOVO_START_IN_FRAME(1),
        /** New start codon is created by the given variant using the chosen transcript. However, it is out of frame relative to the coded protein. */
        DE_NOVO_START_OUT_FRAME(0),

        // These are special / catch-all cases:
        /** Variant lies on one of the RNA transcripts. */
        RNA(4),
        /** Variant lies on one of the lincRNAs. */
        LINCRNA(4);

        /**
         * The relative severity of each {@link VariantClassification}.
         * Lower numbers are considered more severe.
         * Higher numbers are considered less severe.
         */
        final private int relativeSeverity;

        VariantClassification(final int sev) {
            relativeSeverity = sev;
        }
    }
}
