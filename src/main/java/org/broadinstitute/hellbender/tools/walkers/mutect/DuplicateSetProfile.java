package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.UMIReadFilter;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.List;

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
    private static final int INITIAL_LIST_SIZE = 30;
    private final List<GATKRead> currentFamily = new ArrayList<>(INITIAL_LIST_SIZE);
    private UMI currentUmi = null;

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        if (currentUmi == null){
            currentUmi = new UMI(read);
        } else if (! currentUmi.equalsReadUMI(read)){
            getFamilyMetrics();
            currentFamily.clear();
            currentUmi = new UMI(read);
        } else if (currentUmi.equalsReadUMI(read)){
            currentFamily.add(read);
        }

        // TODO: PROBLEM: different molecules can have the same UMI by the probability (1/4^3)^2
        // that is, UMI collision is real, and could happen. Must address that.
        // But for now --- why not just implement the realigning consensus caller?
    }

    private void getFamilyMetrics() {
        for (GATKRead read : currentFamily){
            int d = 3;
        }
    }

    private void something(final GATKRead read, final ReferenceContext referenceContext){
        final byte[] readBases = read.getBases();
        final byte[] referenceBases = referenceContext.getBases();
        final boolean hasIndels = readBases.length != referenceBases.length;
        // OK Got it, ASM walks the ref bases (i.e. skips insertions in reads.)
        final AlignmentStateMachine asm = new AlignmentStateMachine(read);
        asm.stepForwardOnGenome();
        int numMismatches = 0;
        while (!asm.isRightEdge()){
            final CigarElement cigar = asm.getCurrentCigarElement();
            if (cigar.getOperator() == CigarOperator.MATCH_OR_MISMATCH){
                final byte readBase = readBases[asm.getReadOffset()];
                final byte refBase = referenceBases[asm.getGenomeOffset()];
                final int q = 3;
                if (readBase != refBase){
                    numMismatches++;
                }
            } else if (cigar.getOperator() == CigarOperator.DELETION){
                int a = 3;
            }
            asm.stepForwardOnGenome();
        }

        int d = 3;
        asm.stepForwardOnGenome();
        int s = 3;
        while (!asm.isRightEdge()){
            asm.stepForwardOnGenome();
        }

        // ReadUtils.countMismatches(read, getHeaderForReads(), 10, referenceContext.getBases(), 3);
        int e = 3;
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }
}
