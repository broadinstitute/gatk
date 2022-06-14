package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalTextOutputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Count and print to standard output (and optionally to a file) the total number of reads in a SAM/BAM/CRAM file.
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Example</h3>
 *
 * <pre>
 *   gatk CountReads \
 *     -I input_reads.bam
 * </pre>
 */
@CommandLineProgramProperties(
	summary = "Count and print to standard output (and optionally to a file) the total number of reads in a SAM/BAM/CRAM file",
	oneLineSummary = "Count reads in a SAM/BAM/CRAM file",
    programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
@WorkflowProperties
public final class CountReads extends ReadWalker {

    private long count = 0;

    @ArgumentCollection
    final public OptionalTextOutputArgumentCollection out = new OptionalTextOutputArgumentCollection();

    @Override
    public void apply( final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        ++count;
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("CountReads counted " + count + " total reads");
        out.print(count);

        return count;
    }
}
