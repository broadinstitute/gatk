package org.broadinstitute.hellbender.utils.codecs.sampileup;

import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.Feature;

import java.util.List;

/**
 * A tribble feature representing a SAM pileup.
 *
 * Allows intake of both simple (6-column) or extended/consensus (10/13-column) pileups. Simple pileup features will
 * contain only basic information, no observed alleles or variant/genotype inferences, and so shouldn't be used as
 * input for analysis that requires that information.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMPileupFeature implements Feature {
    public enum VariantType { NONE, SNP, INSERTION, DELETION, INDEL }; 

    private String contig;            // genomic location of this genotyped site
    private int start;
    private int stop;

    private char refBaseChar; // what we have set for the reference base (is set to a '*' for indel!)
    private String refBases;        // the reference base sequence according to NCBI; single base for point mutations, deleted bases for  deletions, empty string for insertions

    private String pileupQuals;     // the read base qualities
    private String pileupBases;     // the read bases themselves

    private List<String> observedAlleles = null;    // The sequences of the observed alleles (e.g. {"A","C"} for point mutation or {"","+CC"} for het. insertion
    private VariantType varType = VariantType.NONE;
    private int nNonref = 0; // number of non-reference alleles observed
    private int eventLength = 0; // number of inserted or deleted bases    

    private double consensusScore = 0;
    private double variantScore = 0;

    /**
     * create the pileup feature.  Default protection so that only other classes in this package can create it.
     */
    SAMPileupFeature() {}

    public String getChr() {
        return getContig();
    }

    protected void setChr(String chr) {
        this.contig = chr;
    }

    @Override
    public String getContig() {
        return contig;
    }

    public int getStart() {
        return start;
    }

    protected void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return stop;
    }

    protected void setEnd(int end) {
        this.stop = end;
    }

    public String getQualsAsString()        { return pileupQuals; }

    protected void setPileupQuals(String pileupQuals) {
        this.pileupQuals = pileupQuals;
    }

    /** Returns reference base for point genotypes or '*' for indel genotypes, as a char.
     *
     */
    public char getRef()            { return refBaseChar; }

    protected void setRef(char ref) {
        this.refBaseChar = ref;
    }

    public int size()               { return pileupQuals.length(); }

    /** Returns pile of observed bases over the current genomic location.
     *
     */
    public String getBasesAsString()        { return pileupBases; }

    protected void setPileupBases(String pileupBases) {
        this.pileupBases = pileupBases;
    }

    /** Returns formatted pileup string for the current genomic location as
     * "location: reference_base observed_base_pile observed_qual_pile"
     */
    public String getPileupString()
    {
        if(start == stop)
            return String.format("%s:%d: %s %s %s", getContig(), getStart(), getRef(), getBasesAsString(), getQualsAsString());
        else
            return String.format("%s:%d-%d: %s %s %s", getContig(), getStart(), getEnd(), getRef(), getBasesAsString(), getQualsAsString());
    }

    /**
     * Gets the bases in byte array form.
     * @return byte array of the available bases.
     */
    public byte[] getBases() {
        return StringUtil.stringToBytes(getBasesAsString());
    }

    /**
     * Gets the Phred base qualities without ASCII offset.
     * @return Phred base qualities.
     */
    public byte[] getQuals() {
        byte[] quals = StringUtil.stringToBytes(getQualsAsString());
        for(int i = 0; i < quals.length; i++) quals[i] -= 33;
        return quals;
    }

    /** Returns bases in the reference allele as a String. For point genotypes, the string consists of a single
     * character (reference base). For indel genotypes, the string is empty for insertions into
     * the reference, or consists of deleted bases for deletions.
     *
     * @return reference allele, forward strand
     */
    public String getFWDRefBases() {
        return refBases;
    }

    protected void setRefBases(String refBases) {
        this.refBases = refBases;
    }

    public List<String> getFWDAlleles()  {
        return observedAlleles;
    }

    protected void setFWDAlleles(List<String> alleles) {
        this.observedAlleles = alleles;
    }

    // ----------------------------------------------------------------------
    //
    // What kind of variant are we?
    //
    // ----------------------------------------------------------------------
    public boolean isSNP() { return varType == VariantType.SNP; }
    public boolean isInsertion() { return varType == VariantType.INSERTION; }
    public boolean isDeletion() { return varType == VariantType.DELETION ; }
    public boolean isIndel() { return isInsertion() || isDeletion() || varType == VariantType.INDEL; }
    public boolean isReference()  { return varType == VariantType.NONE; }

    protected void setVariantType(VariantType variantType) {
        this.varType = variantType;
    }

    public double getVariantConfidence() {
        return variantScore;
    }

    protected void setVariantConfidence(double variantScore) {
        this.variantScore = variantScore;
    }

    public boolean isBiallelic() {
        return nNonref  < 2;
    }

    protected void setNumNonRef(int nNonref) {
        this.nNonref = nNonref;
    }

    public double getConsensusConfidence() {
        return consensusScore;
    }

    protected void setConsensusConfidence(double consensusScore) {
        this.consensusScore = consensusScore;
    }

    public int length() {
        return eventLength;
    }

    protected void setLength(int eventLength) {
        this.eventLength = eventLength;
    }

    public boolean isIndelGenotype() {
        return refBaseChar == '*';
    }


    public boolean isPointGenotype() {
        return ! isIndelGenotype();
    }

    /** Implements method required by GenotypeList interface. If this object represents
     * an indel genotype, then it returns itself through this method. If this object is a
     * point genotype, this method returns null.
     * @return
     */
    public SAMPileupFeature getIndelGenotype() {
        if ( isIndelGenotype() ) return this;
        else return null;
    }

    /** Implements method required by GenotypeList interface. If this object represents
     * a point genotype, then it returns itself through this method. If this object is an
     * indel genotype, this method returns null.
     * @return
     */
    public SAMPileupFeature getPointGenotype() {
        if ( isPointGenotype() ) return this;
        else return null;
    }

    /** Returns true if this object \em is an indel genotype (and thus
     * indel genotype is what it only has).
     * @return
     */
    public boolean hasIndelGenotype() {
        return isIndelGenotype();
    }

    /** Returns true if this object \em is a point genotype (and thus
     * point genotype is what it only has.
     * @return
     */
    public boolean hasPointGenotype() {
        return isPointGenotype();
    }

}
