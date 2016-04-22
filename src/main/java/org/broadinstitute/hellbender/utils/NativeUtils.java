package org.broadinstitute.hellbender.utils;

import org.apache.commons.lang3.SystemUtils;

/**
 * Utilities to provide architecture-dependent native functions
 */
public final class NativeUtils {

    private NativeUtils(){}

    /**
     * boolean : returns whether the current platform is PowerPC, LE only
     */
    public static boolean isPPCArchitecture() {
	return SystemUtils.OS_ARCH.equals("ppc64le");
    }
}
