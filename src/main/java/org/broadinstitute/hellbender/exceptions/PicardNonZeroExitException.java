package org.broadinstitute.hellbender.exceptions;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Exception used to propagate non-zero return values from Picard tools.
 */
public class PicardNonZeroExitException extends RuntimeException {
    private static final long serialVersionUID = 0L;
    private final int picardToolReturnCode;

    /**
     *
     * @param picardToolReturnCode the non-zero return code returned from the Picard tool
     */
    public PicardNonZeroExitException(final int picardToolReturnCode) {
        this.picardToolReturnCode = picardToolReturnCode;
        Utils.validateArg(picardToolReturnCode != 0, "Return code must be non-zero");
    }

    /**
     * @return the actual code value returned by the Picard tool
     */
    public int getToolReturnCode() { return picardToolReturnCode; }
}
