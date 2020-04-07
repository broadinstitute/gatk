package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.DuplicateSetWalker;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.consensus.DuplicateSet;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.List;
import java.util.Random;

@CommandLineProgramProperties(
        summary = "Discard a set fraction of duplicate sets from a UMI-grouped ba",
        oneLineSummary = "Discard a set fraction of duplicate sets from a UMI-grouped bam",
        programGroup = ReadDataManipulationProgramGroup.class
)
/**
 * Given a bam grouped by the same unique molecular identifier (UMI), this tool drops a fraction of duplicate sets and returns a new bam.
 * A duplicate set refers to a group of reads whose fragments start at and end at the same coordinate and share the same UMI.
 *
 * The input bam must have been sorted by UMI using FGBio GroupReadsByUmi (http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html).
 *
 * Use this tool to create, for instance, an insilico mixture of duplex-sequenced samples to simulate tumor subclone.
 * Suppose you wish to simulate a tumor sample in which 5% cells share a common set of somatic mutations
 * in addition to ones common to the entire cell population.
 *
 * If you randomly drop 5% of reads in sample A and 95% of reads in sample B and merge the reduced bams,
 * the resulting mixture skews the family-size distribution to the left. Here the family size refers to the
 * number of sequenced duplicate reads that share the same UMI.
 *
 * To see this, take a cancer sample, in which 5% of cells (i.e. a subclone) share a unique set of somatic mutations,
 * that was processed with duplex-UMIs (i.e. UMIs on both adapters) and high rounds of PCR. Suppose we have the sequence-ready
 * libraries of this sample attached to and amplified on the flowcell. Now, sort the flowcell lawn such that the
 * 5% subclone moves near the top of the flowcell. This subclone must have the same family-size distribution as
 * the rest of the flowcell, at about 5% of the library complexity compared to the entire flowcell.
 *
 * Now imagine replacing this subclone with 5% of the *entire* flowcell from another sample prepared and sequenced similarly.
 * The library complexity of these "graft" reads is higher than that of the original, and, consequently, with other parameters
 * such as the number of PCR cycles and sequencing depth fixed, its family distribution would be skewed left---that is, the family size
 * would be smaller than it should be.
 *
 * This tool address the above problem by dropping a set fraction of _duplicate sets_, rather than reads, at random.
 * Implicit in this approach is that a read and its mate are dropped or retained together.
 * While trivial when the input bam is sorted by UMI and query name, this is far from trivial when one attempts
 * to downsample reads naively with a tool like {@link PrintReads}.
 *
 **/
public class DownsampleByDuplicateSet extends DuplicateSetWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "")
    public File outputBam;

    @Argument(fullName = "DS", doc = "This fraction of duplicate sets in the input bam will be retained")
    public double downsamplingRate;

    @Argument(fullName = "keep-duplex-only", doc = "Discard all duplicate sets that don't have duplex evidence")
    public boolean duplexOnly = false;

    private static final int RANDOM_SEED = 142;
    private RandomGenerator rng;
    private static int numFragments;
    private static int numReads;
    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
        outputWriter = createSAMWriter(IOUtils.getPath(outputBam.getAbsolutePath()), false);
    }

    @Override
    public void apply(DuplicateSet duplicateSet, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (filterDuplicateSet(duplicateSet)){
            return;
        }
        if (rng.nextDouble() < downsamplingRate){
            duplicateSet.getReads().forEach(r -> outputWriter.addRead(r));
            numReads += duplicateSet.getReads().size();
            numFragments += 1;
        }
    }

    @Override
    public Object onTraversalSuccess(){
        outputWriter.close();
        logger.info(String.format("Wrote %d reads", numReads));
        logger.info(String.format("Wrote %d fragments", numFragments));
        return "SUCCESS";
    }

    private boolean filterDuplicateSet(final DuplicateSet duplicateSet){
        if (duplicateSet.getReads().size() % 2 == 1){
            // We only keep reads with mates by default, as that's what fgbio GroupByUMI requires.
            logger.info("Duplicate set that contains an unpaired read discarded: " + duplicateSet.getReads().get(0));
            return true;
        }

        // Experiment: only keep duplex
        if (duplexOnly){
            return duplicateSet.isSimplex();
        }

        return false;
    }
}
