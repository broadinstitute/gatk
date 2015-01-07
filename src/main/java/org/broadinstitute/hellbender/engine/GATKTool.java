package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;

/**
 * Base class for all GATK tools. Tool authors that wish to write a "GATK" tool but not use one of
 * the pre-packaged Walker traversals should feel free to extend this class directly. All other
 * GATK tools should extend one of the Walker classes instead.
 */
@CommandLineProgramProperties(usage = "Generic GATK tool", usageShort = "Generic GATK tool", omitFromCommandLine = true)
public abstract class GATKTool extends CommandLineProgram {

    /**
     * Operations performed just prior to the start of traversal. Should be overridden by tool authors
     * who need to process arguments local to their tool or perform other kinds of local initialization.
     *
     * Default implementation does nothing.
     */
    public void onTraversalStart() {}

    /**
     * A complete traversal from start to finish. Tool authors who wish to "roll their own" traversal
     * from scratch can extend this class directly and implement this method. Walker authors should
     * instead extend a Walker class and implement the Walker-appropriate apply() method, since the
     * Walker base classes implement the various kinds of traversals for you.
     */
    public abstract void traverse();

    /**
     * Operations performed immediately after traversal. Should be overridden by tool authors who
     * need to close local resources, etc., after traversal.
     *
     * Default implementation does nothing.
     */
    public void onTraversalDone() {}

    @Override
    protected Object doWork() {
        onTraversalStart();
        traverse();
        onTraversalDone();
        return 0;
    }
}
