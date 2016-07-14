package org.broadinstitute.hellbender.cmdline;

/**
 * only for testing
 */
public final class TestProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return "Testing";
    }

    @Override
    public String getDescription() {
        return "group used for testing";
    }
}