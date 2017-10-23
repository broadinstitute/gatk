package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;

public class BasicValidationResult implements Locatable {
    private int minValidationReadCount;
    private boolean isEnoughValidationReads;
    private boolean isOutOfNoiseFloor;
    private double power;
    private int validationAltCount;
    private int validationRefCount;
    private int discoveryAltCount;
    private int discoveryRefCount;
    private Locatable interval;
    private Allele reference;
    private Allele alternate;
    private String filters;

    public BasicValidationResult(final Locatable interval, int minValidationReadCount, boolean isEnoughValidationReads,
                                 boolean isOutOfNoiseFloor, double power, int validationAltCount,
                                 int validationRefCount, int discoveryAltCount, int discoveryRefCount, final Allele ref,
                                 final Allele alt, final String filters) {
        this.minValidationReadCount = minValidationReadCount;
        this.isEnoughValidationReads = isEnoughValidationReads;
        this.isOutOfNoiseFloor = isOutOfNoiseFloor;
        this.power = power;
        this.validationAltCount = validationAltCount;
        this.validationRefCount = validationRefCount;
        this.discoveryAltCount = discoveryAltCount;
        this.discoveryRefCount = discoveryRefCount;
        this.interval = interval;
        this.alternate = alt;
        this.reference = ref;
        this.filters = filters;
    }

    public int getValidationAltCount() {
        return validationAltCount;
    }

    public int getValidationRefCount() {
        return validationRefCount;
    }

    public int getDiscoveryAltCount() {
        return discoveryAltCount;
    }

    public int getDiscoveryRefCount() {
        return discoveryRefCount;
    }

    public int getMinValidationReadCount() {

        return minValidationReadCount;
    }

    public void setMinValidationReadCount(int minValidationReadCount) {
        this.minValidationReadCount = minValidationReadCount;
    }

    public boolean isEnoughValidationReads() {
        return isEnoughValidationReads;
    }

    public void setEnoughValidationReads(boolean enoughValidationReads) {
        isEnoughValidationReads = enoughValidationReads;
    }

    public boolean isOutOfNoiseFloor() {
        return isOutOfNoiseFloor;
    }

    public void setOutOfNoiseFloor(boolean outOfNoiseFloor) {
        isOutOfNoiseFloor = outOfNoiseFloor;
    }

    public double getPower() {
        return power;
    }

    public void setPower(double power) {
        this.power = power;
    }

    public Locatable getInterval() {
        return interval;
    }

    public void setInterval(Locatable interval) {
        this.interval = interval;
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    public Allele getReference() {
        return reference;
    }

    public void setReference(Allele reference) {
        this.reference = reference;
    }

    public Allele getAlternate() {
        return alternate;
    }

    public void setAlternate(Allele alternate) {
        this.alternate = alternate;
    }

    public String getFilters() {
        return filters;
    }

    public void setFilters(String filters) {
        this.filters = filters;
    }
}
