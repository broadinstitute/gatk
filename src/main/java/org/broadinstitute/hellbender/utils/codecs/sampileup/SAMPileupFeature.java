package org.broadinstitute.hellbender.utils.codecs.sampileup;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;

import java.util.ArrayList;
import java.util.List;

/**
 * A tribble feature representing a SAM pileup.
 *
 * Allows intake of simple mpileups. Simple pileup features will contain only basic information, no reconstructed reads.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMPileupFeature implements Feature {

    private String contig;            // genomic location of this genotyped site
    private int start;

    private char refBaseChar; // what we have set for the reference base (is set to a '*' for indel!)

    private String pileupQuals;     // the read base qualities
    private String pileupBases;     // the read bases themselves

    // TODO: change for PileupElement list
    // private final List<PileupElement> elementList = new ArrayList<>();

    /**
     * create the pileup feature.  Default protection so that only other classes in this package can create it.
     */
    SAMPileupFeature() {}

    @Override
    @Deprecated
    public String getChr() {
        return getContig();
    }

    @Override
    public String getContig() {
        return contig;
    }

    protected void setChr(final String chr) {
        this.contig = chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    protected void setStart(final int start) {
        this.start = start;
    }

    @Override
    public int getEnd() {
        return start;
    }

    public String getQualsAsString() {
        return pileupQuals;
    }

    protected void setPileupQuals(final String pileupQuals) {
        this.pileupQuals = pileupQuals;
    }

    /**
     * Returns reference base for point genotypes or '*' for indel genotypes, as a char.
     */
    public char getRef() {
        return refBaseChar;
    }

    protected void setRef(char ref) {
        this.refBaseChar = ref;
    }

    public int size() {
        return pileupQuals.length();
    }

    /**
     * Returns pile of observed bases over the current genomic location.
     */
    public String getBasesAsString() {
        return pileupBases;
    }

    protected void setPileupBases(String pileupBases) {
        this.pileupBases = pileupBases;
    }

    /**
     * Returns formatted pileup string for the current genomic location as "location: reference_base observed_base_pile
     * observed_qual_pile"
     */
    public String getPileupString() {
        return String.format("%s:%d: %s %s %s",
				getContig(), getStart(), getRef(), getBasesAsString(), getQualsAsString());
    }

    /**
     * Gets the bases in byte array form.
     *
     * @return byte array of the available bases.
     */
    public byte[] getBases() {
        return StringUtil.stringToBytes(getBasesAsString());
    }

    /**
     * Gets the Phred base qualities without ASCII offset.
     *
     * @return Phred base qualities.
     */
    public byte[] getQuals() {
        return SAMUtils.fastqToPhred(getQualsAsString());
    }

}
