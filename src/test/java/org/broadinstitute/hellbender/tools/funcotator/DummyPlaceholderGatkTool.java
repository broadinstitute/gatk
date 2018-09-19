package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;

/**
 * A Dummy / Placeholder class that can be used where a {@link GATKTool} is required.
 * Created by jonn on 9/19/18.
 */
@CommandLineProgramProperties(
        summary = "A dummy GATKTool to help test Funcotator.",
        oneLineSummary = "Dummy dumb dumb tool for testing.",
        programGroup = TestProgramGroup.class
)
public final class DummyPlaceholderGatkTool extends GATKTool {

    @Override
    public void traverse() {

    }

    /**
     * Initialize this {@link DummyPlaceholderGatkTool} by running {@link GATKTool#onStartup()}.
     * @return {@code this} {@link DummyPlaceholderGatkTool}.
     */
    public DummyPlaceholderGatkTool initialize() {
        onStartup();
        return this;
    }
}
