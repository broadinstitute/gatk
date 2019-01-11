package org.broadinstitute.hellbender.cmdline.GATKPlugin.testpluggables;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class TestReadFilter extends ReadFilter {
    private static final long serialVersionUID = 0L;
    @Override
    public boolean test(GATKRead read) {
        return true;
    }
}
