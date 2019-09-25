package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
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
 *   gatk DuplicateSetProfile \
 *     -I input_reads.bam
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = CoverageAnalysisProgramGroup.class
)
public final class DuplicateSetProfile extends ReadWalker {

    private long count = 0;
    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        final AlignmentStateMachine asm = new AlignmentStateMachine(read);
        asm.stepForwardOnGenome();
        for (int i = 0; i < read.getFragmentLength(); i++){
            how exactly do i iterate the read...hmm start here tomorrow 9/23/19 - Mon
        }
    }

    @Override
    public Object onTraversalSuccess() {
        return count;
    }
}
