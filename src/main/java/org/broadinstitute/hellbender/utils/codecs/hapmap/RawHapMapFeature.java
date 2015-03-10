package org.broadinstitute.hellbender.utils.codecs.hapmap;

import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;

import java.util.HashMap;
import java.util.Map;

/**
 * a feature returned by the HapMap Codec - it represents contig, position, name,
 * alleles, other hapmap information, and genotypes for specified samples
 */
public class RawHapMapFeature implements Feature {

    public static final String NULL_ALLELE_STRING = "-";
    public static final String INSERTION = "I";
    public static final String DELETION = "D";

    // the variables we store internally in the class
    private final String name;
    private final String[] alleles;
    private Map<String, Allele> actualAlleles = null;
    private final String contig ;
    private long position;
    private final Strand strand;
    private final String assembly;
    private final String center;
    private final String protLSID;
    private final String assayLSID;
    private final String panelLSID;
    private final String qccode;
    private final String[] genotypes;

    // we store the header line, if they'd like to get the samples
    private final String headerLine;

    /**
     * create a HapMap Feature, based on all the records available in the hapmap file
     * @param contig the contig name
     * @param position the position
     * @param strand the strand enum
     * @param assembly what assembly this feature is from
     * @param center the center that provided this SNP
     * @param protLSID ??
     * @param assayLSID ??
     * @param panelLSID ??
     * @param qccode ??
     * @param genotypes a list of strings, representing the genotypes for the list of samples
     */
    public RawHapMapFeature(String name,
                            String[] alleles,
                            String contig,
                            Long position,
                            Strand strand,
                            String assembly,
                            String center,
                            String protLSID,
                            String assayLSID,
                            String panelLSID,
                            String qccode,
                            String[] genotypes,
                            String headerLine) {
        this.name = name;
        this.alleles = alleles;
        this.contig = contig;
        this.position = position;
        this.strand = strand;
        this.assembly =  assembly;
        this.center =  center;
        this.protLSID = protLSID ;
        this.assayLSID = assayLSID ;
        this.panelLSID = panelLSID ;
        this.qccode = qccode;
        this.genotypes = genotypes;
        this.headerLine = headerLine;
    }

    /**
     * get the contig value
     * @return a string representing the contig
     */
    public String getChr() {
        return contig;
    }

    /**
     * get the start position, as an integer
     * @return an int, representing the start position
     */
    public int getStart() {
        return (int)position;
    }

    /**
     * get the end position
     * @return get the end position as an int
     */
    public int getEnd() {
        return (int)position;
    }

    /**
     * Getter methods
     */

    public Strand getStrand() {
        return strand;
    }

    public String getAssembly() {
        return assembly;
    }

    public String getCenter() {
        return center;
    }

    public String getProtLSID() {
        return protLSID;
    }

    public String getAssayLSID() {
        return assayLSID;
    }

    public String getPanelLSID() {
        return panelLSID;
    }

    public String getQCCode() {
        return qccode;
    }

    public String getName() {
        return name;
    }

    public String[] getAlleles() {
        return alleles;
    }

    public String[] getGenotypes() {
        return genotypes;
    }

    // This is necessary because HapMap places insertions in the incorrect position
    public void updatePosition(long position) {
        this.position = position;
    }

    public void setActualAlleles(Map<String, Allele> alleleMap) {
        actualAlleles = new HashMap<String, Allele>(alleleMap);
    }

    public Map<String, Allele> getActualAlleles() {
        return actualAlleles;
    }
    
    /**
     * get a list of the samples from the header (in order)
     * @return a string array of sample names
     */
    public String[] getSampleIDs() {
		String[] header = headerLine.split("\\s+");
		String[] sample_ids = new String[header.length-11];
		for (int i = 11; i < header.length; i++)
			sample_ids[i-11] = header[i];
		return sample_ids;
	}
}
