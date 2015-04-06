package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * Illumina's TileMetricsOut.bin file codes various metrics, both concrete (all density id's are code 100) or as a base code
 * (e.g. phasing values are computed from a base of 200).
 *
 * @author jgentry
 */
public enum IlluminaMetricsCode {
    DENSITY_ID(100),
    CLUSTER_ID(102),
    PHASING_BASE(200),
    PREPHASING_BASE(201);

    private final int metricsCode;

    IlluminaMetricsCode(final int metricsCode) {
        this.metricsCode = metricsCode;
    }

    /**
     * Phasing codes are between 200 and 299 (inclusive). Phasing codes are defined as being
     * (200 + ((N - 1) * 2)) for (a 0-based) read descriptor N (i.e., 200, 202, 204, etc.) Prephasing codes are defined
     * as being (201 + ((N - 1) * 2)) for read descriptor N (i.e., 201, 203, 205, etc.). So for a 101T8B101T read
     * structure, there will be phasing codes of 200, 202 and 204 and prephasing codes of 201, 203, 205.
     */
    public static int getPhasingCode(final int readDescriptorIndex, final IlluminaMetricsCode phasingType) {
        if (!isPhasing(phasingType)) {
            throw new IllegalArgumentException("phasingType must be PHASING_BASE or PREPHASING_BASE");
        }

        return (phasingType.getMetricsCode() + (readDescriptorIndex * 2));
    }

    public static boolean isPhasing(final IlluminaMetricsCode metricsCode) {
        return (metricsCode.equals(PHASING_BASE) || metricsCode.equals(PREPHASING_BASE));
    }

    public int getMetricsCode() {
        return metricsCode;
    }
}
