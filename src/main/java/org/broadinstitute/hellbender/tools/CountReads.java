package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Optional;

@CommandLineProgramProperties(
	usage = "Count reads.",
	usageShort = "Count reads",
    programGroup = ReadProgramGroup.class
)
public class CountReads extends ReadWalker {

    private long count = 0;
    @Override
    public void apply( SAMRecord read, Optional<ReferenceContext> referenceContext, Optional<FeatureContext> featureContext ) {
        ++count;
    }

    @Override
    public Object onTraversalDone() {
        return count;
    }
}
