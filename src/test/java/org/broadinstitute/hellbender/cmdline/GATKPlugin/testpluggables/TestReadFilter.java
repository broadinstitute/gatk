package org.broadinstitute.hellbender.cmdline.GATKPlugin.testpluggables;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * This is a test ReadFilter that is used to test if the gatk config file properly controls ReadFilter loading
 * see {@link org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptorTest}
 *
 * This should not be discoverable normally since it's not in the expected ReadFilter package.
 */
public class TestReadFilter extends ReadFilter {
    private static final long serialVersionUID = 0L;
    @Override
    public boolean test(GATKRead read) {
        return true;
    }
}
