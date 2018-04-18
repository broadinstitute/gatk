package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;

/**
 * Implementation of {@link ShardingArgumentCollection} with optional values for all arguments.
 *
 * <p>Note: this class use the recommended argument names defined in {@link ShardingArgumentCollection}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class OptionalShardingArgumentCollection implements ShardingArgumentCollection {

    @Argument(fullName = WINDOW_SIZE_NAME, doc = "Window-size for the analysis", optional = true)
    public int windowSize;

    @Argument(fullName = WINDOW_STEP_NAME, doc = "Window-step for the analysis", optional = true)
    public int windowStep;

    @Argument(fullName = WINDOW_PAD_NAME, doc = "Window-pad for the analysis", optional = true)
    public int windowPad;

    /**
     * Creates a new argument collection with default values.
     *
     * @param defaultWinSize non-zero positive window-size (default).
     * @param defaultWinStep non-zero positive window-step (default).
     * @param defaultWinPad  positive window-pad (default).
     */
    public OptionalShardingArgumentCollection(final int defaultWinSize, final int defaultWinStep,
            final int defaultWinPad) {
        this.windowSize = defaultWinSize;
        this.windowStep = defaultWinStep;
        this.windowPad = defaultWinPad;
        try {
            validate();
        } catch (final CommandLineException e) {
            throw new IllegalArgumentException(e);
        }
    }

    @Override
    public int getWindowSize() {
        return windowSize;
    }

    @Override
    public int getWindowStep() {
        return windowStep;
    }

    @Override
    public int getWindowPad() {
        return windowPad;
    }
}
