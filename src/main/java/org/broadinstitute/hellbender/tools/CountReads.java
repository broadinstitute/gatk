package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Calculate and print to the standard output the overall number of reads in a SAM/BAM/CRAM file
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
@DocumentedFeature
@CommandLineProgramProperties(
	summary = "Count reads in a SAM/BAM/CRAM file.",
	oneLineSummary = "Count reads in a SAM/BAM/CRAM file",
    programGroup = CoverageAnalysisProgramGroup.class
)
public final class CountReads extends ReadWalker {

    private long count = 0;
    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        ++count;
    }

    @Override
    public Object onTraversalSuccess() {
        return count;
    }
}
