package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;

/**
 * Argument collection for window slider parameters (window-size, window-step, window-padding).
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class SliderArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private WindowSizeArgumentCollection windowSize;

    @ArgumentCollection
    private WindowStepArgumentCollection windowStep;

    @ArgumentCollection
    private WindowPaddingArgumentCollection windowPadding;

    public SliderArgumentCollection(final int defaultWindowSize, final boolean hasWindowSizeArgument,
                                    final int defaultWindowStep, final boolean hasWindowStepArgument,
                                    final int defaultWindowPadding, final boolean hasWindowPaddingArgument) {
        // set the window size
        windowSize = (hasWindowSizeArgument) ?
                new OptionalWindowSizeArgumentCollection(defaultWindowSize) :
                new DefaultWindowSizeArgumentCollection(defaultWindowSize);
        // set the window step
        windowStep = (hasWindowStepArgument) ?
                new OptionalWindowStepArgumentCollection(defaultWindowStep) :
                new DefaultWindowStepArgumentCollection(defaultWindowStep);
        // set the window padding argument
        windowPadding = (hasWindowPaddingArgument) ?
                new OptionalWindowPaddingArgumentCollection(defaultWindowPadding) :
                new DefaultWindowPaddingArgumentCollection(defaultWindowPadding);
    }

    public int getWindowSize() {
        return windowSize.getWindowSize();
    }

    public int getWindowStep() {
        return windowStep.getWindowStep();
    }

    public int getWindowPadding() {
        return windowPadding.getWindowPadding();
    }

}
