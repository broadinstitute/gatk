package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;

public class BasicValidationResult implements Locatable {
    private final int minValidationReadCount;
    private final boolean isEnoughValidationReads;
    private final boolean isOutOfNoiseFloor;
    private final double power;
    private final int validationAltCount;
    private final int validationRefCount;
    private final int discoveryAltCount;
    private final int discoveryRefCount;
    private final Locatable interval;
    private final Allele reference;
    private final Allele alternate;
    private final String filters;
    private final long numAltSupportingReadsInNormal;

    public BasicValidationResult(final Locatable interval, final int minValidationReadCount, final boolean isEnoughValidationReads,
                                 final boolean isOutOfNoiseFloor, final double power, final int validationAltCount,
                                 final int validationRefCount, final int discoveryAltCount, final int discoveryRefCount, final Allele ref,
                                 final Allele alt, final String filters, final long numAltSupportingReadsInNormal) {
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
        this.numAltSupportingReadsInNormal = numAltSupportingReadsInNormal;
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

    public boolean isEnoughValidationReads() {
        return isEnoughValidationReads;
    }

    public boolean isOutOfNoiseFloor() {
        return isOutOfNoiseFloor;
    }

    public double getPower() {
        return power;
    }

    public Locatable getInterval() {
        return interval;
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

    public Allele getAlternate() {
        return alternate;
    }

    public String getFilters() {
        return filters;
    }

    public long getNumAltSupportingReadsInNormal() {
        return numAltSupportingReadsInNormal;
    }
}
