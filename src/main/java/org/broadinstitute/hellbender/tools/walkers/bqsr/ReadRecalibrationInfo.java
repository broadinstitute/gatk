package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

public final class ReadRecalibrationInfo {
    private final SAMRecord read;
    private final int length;
    private final ReadCovariates covariates;
    private final boolean[] skips;
    private final byte[] baseQuals, insertionQuals, deletionQuals;
    private final double[] snpErrors, insertionErrors, deletionErrors;

    public ReadRecalibrationInfo(final SAMRecord read,
                                 final ReadCovariates covariates,
                                 final boolean[] skips,
                                 final double[] snpErrors,
                                 final double[] insertionErrors,
                                 final double[] deletionErrors) {
        if ( read == null ) throw new IllegalArgumentException("read cannot be null");
        if ( covariates == null ) throw new IllegalArgumentException("covariates cannot be null");
        if ( skips == null ) throw new IllegalArgumentException("skips cannot be null");
        if ( snpErrors == null ) throw new IllegalArgumentException("snpErrors cannot be null");
        if ( insertionErrors == null ) throw new IllegalArgumentException("insertionErrors cannot be null");
        if ( deletionErrors == null ) throw new IllegalArgumentException("deletionErrors cannot be null");

        this.read = read;
        this.baseQuals = read.getBaseQualities();
        this.length = baseQuals.length;
        this.covariates = covariates;
        this.skips = skips;
        this.insertionQuals = ReadUtils.getExistingBaseInsertionQualities(read);
        this.deletionQuals = ReadUtils.getExistingBaseDeletionQualities(read);
        this.snpErrors = snpErrors;
        this.insertionErrors = insertionErrors;
        this.deletionErrors = deletionErrors;

        if ( skips.length != length ) throw new IllegalArgumentException("skips.length " + snpErrors.length + " != length " + length);
        if ( snpErrors.length != length ) throw new IllegalArgumentException("snpErrors.length " + snpErrors.length + " != length " + length);
        if ( insertionErrors.length != length ) throw new IllegalArgumentException("insertionErrors.length " + snpErrors.length + " != length " + length);
        if ( deletionErrors.length != length ) throw new IllegalArgumentException("deletionErrors.length " + snpErrors.length + " != length " + length);
    }

    /**
     * Get the qual score for event type at offset
     *
     * @param eventType the type of event we want the qual for
     * @param offset the offset into this read for the qual
     * @return a valid quality score for event at offset
     */
    public byte getQual(final EventType eventType, final int offset) {
        switch ( eventType ) {
            case BASE_SUBSTITUTION: return baseQuals[offset];
            // note optimization here -- if we don't have ins/del quals we just return the default byte directly
            case BASE_INSERTION: return insertionQuals == null ? ReadUtils.DEFAULT_INSERTION_DELETION_QUAL : insertionQuals[offset];
            case BASE_DELETION: return deletionQuals == null ? ReadUtils.DEFAULT_INSERTION_DELETION_QUAL : deletionQuals[offset];
            default: throw new IllegalStateException("Unknown event type " + eventType);
        }
    }

    /**
     * Get the error fraction for event type at offset
     *
     * The error fraction is a value between 0 and 1 that indicates how much certainty we have
     * in the error occurring at offset.  A value of 1 means that the error definitely occurs at this
     * site, a value of 0.0 means it definitely doesn't happen here.  0.5 means that half the weight
     * of the error belongs here
     *
     * @param eventType the type of event we want the qual for
     * @param offset the offset into this read for the qual
     * @return a fractional weight for an error at this offset
     */
    public double getErrorFraction(final EventType eventType, final int offset) {
        switch ( eventType ) {
            case BASE_SUBSTITUTION: return snpErrors[offset];
            case BASE_INSERTION: return insertionErrors[offset];
            case BASE_DELETION: return deletionErrors[offset];
            default: throw new IllegalStateException("Unknown event type " + eventType);
        }
    }

    /**
     * Get the read involved in this recalibration info
     * @return a non-null GATKSAMRecord
     */
    public SAMRecord getRead() {
        return read;
    }

    /**
     * Should offset in this read be skipped (because it's covered by a known variation site?)
     * @param offset a valid offset into this info
     * @return true if offset should be skipped, false otherwise
     */
    public boolean skip(final int offset) {
        return skips[offset];
    }

    /**
     * Get the ReadCovariates object carrying the mapping from offsets -> covariate key sets
     * @return a non-null ReadCovariates object
     */
    public ReadCovariates getCovariatesValues() {
        return covariates;
    }

    /**
     * Ensures an offset is valid.  Used in contracts
     * @param offset a proposed offset
     * @return true if offset is valid w.r.t. the data in this object, false otherwise
     */
    private boolean validOffset(final int offset) {
        return offset >= 0 && offset < baseQuals.length;
    }

    private boolean validQual(final byte result) {
        return result >= 0 && result <= QualityUtils.MAX_SAM_QUAL_SCORE;
    }
}
