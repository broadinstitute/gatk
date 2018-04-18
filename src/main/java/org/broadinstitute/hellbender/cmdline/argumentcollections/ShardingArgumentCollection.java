package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.CommandLineException;

/**
 * Interface for arguments related with sharding intervals over windows.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public interface ShardingArgumentCollection {

    /** Recommended valur for window-size argument. */
    public static final String WINDOW_SIZE_NAME = "window-size";
    /** Recommended valur for window-step argument. */
    public static final String WINDOW_STEP_NAME = "window-step";
    /** Recommended valur for window-pad argument. */
    public static final String WINDOW_PAD_NAME = "window-pad";

    /**
     * Returns the window-size applied to the data.
     */
    public int getWindowSize();

    /**
     * Returns the window-step applied to the data.
     *
     * <p>If equal to {@link #getWindowSize()}, the analysis would be performed in non-overlapping windows.
     */
    public int getWindowStep();

    /**
     * Returns the window-padding applied to the data.
     *
     * <p>If it returns {@code 0}, no padding is expected.
     */
    public int getWindowPad();

    /**
     * Validates arguments.
     *
     * <p>Default method checks that window-size and window-step are larger than 1 and
     * that window-pad is not negative.
     *
     * @throws CommandLineException if validation fails.
     */
    public default void validate() {
        if (getWindowSize() < 1) {
            throw new CommandLineException.BadArgumentValue(WINDOW_SIZE_NAME,
                    String.valueOf(getWindowSize()), "should be positive (non-zero)");
        }
        if (getWindowStep() < 1) {
            throw new CommandLineException.BadArgumentValue(WINDOW_STEP_NAME,
                    String.valueOf(getWindowStep()), "should be positive (non-zero)");
        }
        if (getWindowPad() < 0) {
            throw new CommandLineException.BadArgumentValue(WINDOW_PAD_NAME,
                    String.valueOf(getWindowPad()), "should be positive (or zero)");
        }
    }
}
