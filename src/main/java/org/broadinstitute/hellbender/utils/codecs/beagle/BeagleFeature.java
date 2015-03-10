package org.broadinstitute.hellbender.utils.codecs.beagle;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;

import java.util.ArrayList;
import java.util.Map;

public class BeagleFeature implements Feature {

    private String chr;
    private int start;
    private int end;

    Map<String, ArrayList<String>> sampleGenotypes;
    private Double r2Value;
    Map<String, ArrayList<String>> probLikelihoods;

    Allele AlleleA;
    Allele AlleleB;


    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public Double getR2value() {
        return r2Value;
    }

    public Allele getAlleleA() {
        return AlleleA;
    }

    public Allele getAlleleB() {
        return AlleleB;
    }

    public Map<String, ArrayList<String>> getProbLikelihoods() {
        return probLikelihoods;
    }

    public Map<String, ArrayList<String>> getGenotypes() {
        return sampleGenotypes;        
    }

    protected void setChr(String chr) {
       this.chr = chr;
    }

    protected void setStart(int start) {
        this.start = start;
    }

    protected void setEnd(int end) {
        this.end = end;
    }

    protected void setR2value(double r2) {
        this.r2Value = r2;
    }

    protected void setAlleleA(String a, boolean isRef) {
        this.AlleleA = Allele.create(a, isRef);
    }

    protected void setAlleleB(String a, boolean isRef) {
        this.AlleleB = Allele.create(a, isRef);
    }

    protected void setProbLikelihoods(Map<String, ArrayList<String>> p) {
        this.probLikelihoods = p;
    }

    protected void setGenotypes(Map<String, ArrayList<String>> p) {
        this.sampleGenotypes = p;
    }
}
