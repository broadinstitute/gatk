package org.broadinstitute.hellbender.tools.copynumber.formats;

public final class CopyNumberFormatsUtils {
    public static final String DOUBLE_FORMAT = "%.6f";

    private CopyNumberFormatsUtils() {}

    public static String formatDouble(final double value) {
        return String.format(DOUBLE_FORMAT, value);
    }
}
