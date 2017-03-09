package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
// TODO: maybe this should be merged with ReadFilterPluginUnitTest
public class GATKReadFilterPluginDescriptorUnitTest extends BaseTest {

    @CommandLineProgramProperties(summary = "Test read filter plugin with default arguments",
            oneLineSummary = "Test read filter plugin with default arguments",
            programGroup = TestProgramGroup.class)
    private static class TestWithDefaultReadFilters extends CommandLineProgram {

        private final List<ReadFilter> defaultFilters;

        public TestWithDefaultReadFilters(final List<ReadFilter> defaultFilters) {
            this.defaultFilters = defaultFilters;
        }

        protected List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
            return Collections.singletonList(new GATKReadFilterPluginDescriptor(defaultFilters));
        }

        @Override
        protected Object doWork() {
            return null;
        }
    }

    @Test
    public void testWithDefaultReadFiltersWithParams() throws Exception {
        // this ReadFilter have parameters --maxReadLength/--minReadLength, that are set because of the default filter
        final CommandLineProgram clp = new TestWithDefaultReadFilters(Collections.singletonList(new ReadLengthReadFilter(10, 50)));
        // disable the read filter should not blow up because of that parameters, because they are not provided by the user
        clp.instanceMain(new String[]{"--" + StandardArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME, "ReadLengthReadFilter"});

    }
}