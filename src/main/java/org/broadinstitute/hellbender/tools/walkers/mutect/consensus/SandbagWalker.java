package org.broadinstitute.hellbender.tools.walkers.mutect.consensus;

import htsjdk.samtools.CigarElement;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.kohsuke.args4j.Argument;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

/** Test whatever your heart desires here **/
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "Print reads in the SAM/BAM/CRAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
public class SandbagWalker extends ReadWalker {

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        final AlignmentStateMachine asm = new AlignmentStateMachine(read);
        for (int i = 0; i < read.getLength(); i++){
            asm.stepForwardOnGenome();
            CigarElement ce = asm.getCurrentCigarElement();
        }

    }
}
