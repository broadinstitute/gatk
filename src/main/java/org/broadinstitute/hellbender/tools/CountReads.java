package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

@CommandLineProgramProperties(
	summary = "Count reads.",
	oneLineSummary = "Count reads",
    programGroup = ReadProgramGroup.class
)
public final class CountReads extends ReadWalker {

    private long count = 0;
    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        ++count;
    }

    @Override
    public Object onTraversalDone() {
        return count;
    }
}
