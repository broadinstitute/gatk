package org.broadinstitute.hellbender.tools.gvs.common;

import java.util.Objects;

/**
 * Canonical VRS allele record.
 *
 * <p>One instance is produced per normalized alternate allele during VET ingest.
 * It carries every column needed to populate the {@code vrs_allele} table:
 * <ul>
 *   <li>{@link #vrsAlleleId} — {@code ga4gh:VA.*} (also stored as FK in VET vrs_allele_ids)</li>
 *   <li>{@link #vrsLocationId} — {@code ga4gh:SL.*} SequenceLocation identifier</li>
 *   <li>{@link #refgetAccession} — GA4GH RefGet {@code SQ.*} chromosome accession</li>
 *   <li>{@link #start} / {@link #end} — 0-based interbase coordinates (VRS)</li>
 *   <li>{@link #stateType} — {@code LiteralSequenceExpression} or {@code ReferenceLengthExpression}</li>
 *   <li>{@link #stateSequence} — alt sequence for literal states; null for RLE</li>
 *   <li>{@link #stateLength} — total alt length for RLE; null for literal</li>
 *   <li>{@link #stateRepeatSubunitLength} — repeat unit length for RLE; null for literal</li>
 * </ul>
 */
public class VrsAlleleRecord {

    public static final String STATE_LITERAL   = "LiteralSequenceExpression";
    public static final String STATE_REF_LEN   = "ReferenceLengthExpression";

    /** {@code ga4gh:VA.<digest>} */
    public final String vrsAlleleId;
    /** {@code ga4gh:SL.<digest>} */
    public final String vrsLocationId;
    /** GA4GH RefGet accession, e.g. {@code SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl} */
    public final String refgetAccession;
    /** 0-based interbase start (VRS coordinates) */
    public final long start;
    /** 0-based interbase end (VRS coordinates, half-open) */
    public final long end;
    /** {@code LiteralSequenceExpression} or {@code ReferenceLengthExpression} */
    public final String stateType;
    /** Alternate sequence for literal states; null for ReferenceLengthExpression */
    public final String stateSequence;
    /** Total alt allele length for RLE states; null for LiteralSequenceExpression */
    public final Integer stateLength;
    /** Repeat unit length for RLE states; null for LiteralSequenceExpression */
    public final Integer stateRepeatSubunitLength;

    /** Constructor for a LiteralSequenceExpression allele (SNPs, MNPs, unambiguous insertions). */
    public VrsAlleleRecord(String vrsAlleleId, String vrsLocationId, String refgetAccession,
                           long start, long end, String stateSequence) {
        this.vrsAlleleId              = vrsAlleleId;
        this.vrsLocationId            = vrsLocationId;
        this.refgetAccession          = refgetAccession;
        this.start                    = start;
        this.end                      = end;
        this.stateType                = STATE_LITERAL;
        this.stateSequence            = stateSequence;
        this.stateLength              = null;
        this.stateRepeatSubunitLength = null;
    }

    /** Constructor for a ReferenceLengthExpression allele (deletions, tandem repeats). */
    public VrsAlleleRecord(String vrsAlleleId, String vrsLocationId, String refgetAccession,
                           long start, long end, int stateLength, int stateRepeatSubunitLength) {
        this.vrsAlleleId              = vrsAlleleId;
        this.vrsLocationId            = vrsLocationId;
        this.refgetAccession          = refgetAccession;
        this.start                    = start;
        this.end                      = end;
        this.stateType                = STATE_REF_LEN;
        this.stateSequence            = null;
        this.stateLength              = stateLength;
        this.stateRepeatSubunitLength = stateRepeatSubunitLength;
    }

    @Override
    public String toString() {
        return String.format("VrsAlleleRecord{vrsAlleleId='%s', vrsLocationId='%s', "
                + "refgetAccession='%s', start=%d, end=%d, stateType='%s'}",
                vrsAlleleId, vrsLocationId, refgetAccession, start, end, stateType);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof VrsAlleleRecord)) return false;
        VrsAlleleRecord that = (VrsAlleleRecord) o;
        return start == that.start
            && end == that.end
            && Objects.equals(vrsAlleleId, that.vrsAlleleId)
            && Objects.equals(vrsLocationId, that.vrsLocationId)
            && Objects.equals(refgetAccession, that.refgetAccession)
            && Objects.equals(stateType, that.stateType)
            && Objects.equals(stateSequence, that.stateSequence)
            && Objects.equals(stateLength, that.stateLength)
            && Objects.equals(stateRepeatSubunitLength, that.stateRepeatSubunitLength);
    }

    @Override
    public int hashCode() {
        return Objects.hash(vrsAlleleId, vrsLocationId, refgetAccession, start, end,
                stateType, stateSequence, stateLength, stateRepeatSubunitLength);
    }
}
