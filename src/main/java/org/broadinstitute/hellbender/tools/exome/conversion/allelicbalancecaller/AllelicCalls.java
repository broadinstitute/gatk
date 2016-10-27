package org.broadinstitute.hellbender.tools.exome.conversion.allelicbalancecaller;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;

import java.io.Serializable;

public class AllelicCalls implements Serializable, Locatable {

    static final long serialVersionUID = 337337337L;


    private ACNVModeledSegment acnvSegment;
    private AllelicBalanceCall balancedCall;
    private CNLoHCall cnlohCall;
    private double fCr;
    private double fMaf;
    private double rho;
    private int m;
    private int n;


    public AllelicCalls(final ACNVModeledSegment acnvSegment) {
        this.acnvSegment = acnvSegment;
        balancedCall = AllelicBalanceCall.NO_CALL;
        cnlohCall = CNLoHCall.NO_CALL;
    }

    public ACNVModeledSegment getAcnvSegment() {
        return acnvSegment;
    }

    public void setAcnvSegment(final ACNVModeledSegment acnvSegment) {
        this.acnvSegment = acnvSegment;
    }

    public AllelicBalanceCall getBalancedCall() {
        return balancedCall;
    }

    public void setBalancedCall(final AllelicBalanceCall balancedCall) {
        this.balancedCall = balancedCall;
    }

    public CNLoHCall getCnlohCall() {
        return cnlohCall;
    }

    public void setCnlohCall(final CNLoHCall cnlohCall) {
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

    public void setfCr(final double fCr) {
        this.fCr = fCr;
    }

    public double getfMaf() {
        return fMaf;
    }

    public void setfMaf(final double fMaf) {
        this.fMaf = fMaf;
    }

    public double getRho() {
        return rho;
    }

    public void setRho(final double rho) {
        this.rho = rho;
    }

    public int getM() {
        return m;
    }

    public void setM(final int m) {
        this.m = m;
    }

    public int getN() {
        return n;
    }

    public void setN(final int n) {
        this.n = n;
    }
}
