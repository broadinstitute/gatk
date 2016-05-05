package org.broadinstitute.hellbender.tools.exome.cnlohcaller;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;

import java.io.Serializable;

public class CNLOHCall implements Serializable, Locatable{

    static final long serialVersionUID = 337337337L;


    private ACNVModeledSegment acnvSegment;
    private CNLOHBalancedCallEnum balancedCall;
    private CNLOHLoHCallEnum cnlohCall;
    private double fCr;
    private double fMaf;
    private double rho;
    private int m;
    private int n;


    public CNLOHCall(final ACNVModeledSegment acnvSegment) {
        this.acnvSegment = acnvSegment;
        balancedCall = CNLOHBalancedCallEnum.NO_CALL;
        cnlohCall = CNLOHLoHCallEnum.NO_CALL;
    }


    public ACNVModeledSegment getAcnvSegment() {
        return acnvSegment;
    }

    public void setAcnvSegment(final ACNVModeledSegment acnvSegment) {
        this.acnvSegment = acnvSegment;
    }

    public CNLOHBalancedCallEnum getBalancedCall() {
        return balancedCall;
    }

    public void setBalancedCall(final CNLOHBalancedCallEnum balancedCall) {
        this.balancedCall = balancedCall;
    }

    public CNLOHLoHCallEnum getCnlohCall() {
        return cnlohCall;
    }

    public void setCnlohCall(final CNLOHLoHCallEnum cnlohCall) {
        this.cnlohCall = cnlohCall;
    }

    @Override
    public String getContig() {
        return acnvSegment.getContig();
    }

    @Override
    public int getStart() {
        return acnvSegment.getStart();
    }

    @Override
    public int getEnd() {
        return acnvSegment.getEnd();
    }

    public double getfCr() {
        return fCr;
    }

    public void setfCr(double fCr) {
        this.fCr = fCr;
    }

    public double getfMaf() {
        return fMaf;
    }

    public void setfMaf(double fMaf) {
        this.fMaf = fMaf;
    }

    public double getRho() {
        return rho;
    }

    public void setRho(double rho) {
        this.rho = rho;
    }

    public int getM() {
        return m;
    }

    public void setM(int m) {
        this.m = m;
    }

    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }
}
