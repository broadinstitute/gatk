package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.CommandLineException;

/**
 * Interface for arguments related to sharding intervals over windows.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public interface ShardingArgumentCollection {

    /** Recommended value for window-size argument. */
    public static final String WINDOW_SIZE_NAME = "window-size";
    /** Recommended value for window-step argument. */
    public static final String WINDOW_STEP_NAME = "window-step";
    /** Recommended value for window-pad argument. */
    public static final String WINDOW_PADDING_NAME = "window-padding";

    /**
     * Returns the number of bases included in each windows (window size).
     *
     * @return window size applied to the data.
     */
    public int getWindowSize();

    /**
     * Returns the number of bases overlapping between consecutive windows (window step).
     *
     * <p>If equal to {@link #getWindowSize()}, the analysis would be performed in non-overlapping windows.
     *
     * @return window step applied to the data.
     */
    public int getWindowStep();

    /**
     * Returns the number of bases to be padded on both sides of each window (window padding).
     *
     * <p>If it returns {@code 0}, no padding is expected.
     *
     * @return window padding applied to the data.
     */
    public int getWindowPadding();

    /**
     * Validates arguments.
     *
     * <p>Default method checks that window size and step are larger than 1 and
     * that window padding is not negative.
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
        if (getWindowPadding() < 0) {
            throw new CommandLineException.BadArgumentValue(WINDOW_PADDING_NAME,
                    String.valueOf(getWindowPadding()), "should be positive (or zero)");
        }
    }
}
