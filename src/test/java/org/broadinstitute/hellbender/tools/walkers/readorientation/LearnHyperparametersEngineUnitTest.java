package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.Histogram;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.ArtifactState;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.CollectDataForReadOrientationFilter.MAX_REF_DEPTH;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.CollectDataForReadOrientationFilter.REGULAR_BASES;

/**R
 * Created by tsato on 8/14/17.
 */
public class LearnHyperparametersEngineUnitTest {
    private static final double EPSILON = 1e-2;

    final String contig = "1";
    final int position = 100_000_000;
    final int DEPTH = 100;
    final int ALT_DEPTH = 20;
    final int[] altSiteDepthCounts = new int[]{0,0,0, ALT_DEPTH};
    final int[] biasedAltF1R2DepthCounts = new int[]{0,0,0, ALT_DEPTH};

    protected final Logger logger = LogManager.getLogger(this.getClass());


    /**
     * 100 sites for a single context AGT. G -> T transition observed in NUM_F1R2_EXAMPLES sites
     * All of the alt sites have 100% Alt F1R2 reads. The rest of the examples are hom ref.
     */
    @Test
    public void testSimpleCase() {
        final int NUM_EXAMPLES = 100;
        final int NUM_ALT_EXAMPLES = 20;
        final int NUM_F1R2_EXAMPLES = NUM_ALT_EXAMPLES;

        final String refContext = "AGT";

        final List<AltSiteRecord> altDeisgnMatrix = new ArrayList<>();
        final Histogram<Integer> refSiteHistogram = new Histogram<>("depth", refContext);
        refSiteHistogram.prefillBins(IntStream.rangeClosed(1, MAX_REF_DEPTH).boxed().toArray(Integer[]::new));

        for (int n = 0; n < NUM_EXAMPLES; n++){
            // assume 20 in 100 examples is alt at ALT_DEPTH/DEPTH % allele fraction, and alt reads are entirely F1R2
            if (n < NUM_ALT_EXAMPLES){
                altDeisgnMatrix.add(new AltSiteRecord(contig, position, refContext, altSiteDepthCounts,
                        biasedAltF1R2DepthCounts, DEPTH, Nucleotide.T));
            } else {
                refSiteHistogram.increment(DEPTH);
            }
        }

        LearnHyperparametersEngine engine = new LearnHyperparametersEngine(refSiteHistogram, altDeisgnMatrix,
                LearnHyperparameters.DEFAULT_CONVERGENCE_THRESHOLD, LearnHyperparameters.DEFAULT_MAX_ITERATIONS, logger);
        Hyperparameters hyperparameters = engine.runEMAlgorithm();

        Assert.assertEquals(engine.effectiveCounts[ArtifactState.F1R2_T.ordinal()], (double) NUM_F1R2_EXAMPLES, EPSILON);
        Assert.assertEquals(engine.effectiveCounts[ArtifactState.HOM_REF.ordinal()], (double) NUM_EXAMPLES - NUM_ALT_EXAMPLES, EPSILON);

        // check pi
        Assert.assertEquals(hyperparameters.getPi()[LearnHyperparametersEngine.ArtifactState.F1R2_T.ordinal()], (double) NUM_ALT_EXAMPLES/NUM_EXAMPLES, EPSILON);
        Assert.assertEquals(hyperparameters.getPi()[LearnHyperparametersEngine.ArtifactState.HOM_REF.ordinal()], 1.0 - (double) NUM_ALT_EXAMPLES/NUM_EXAMPLES, EPSILON);
    }


    /**
     * Now test the case where not all of the transitions have orientation bias. And for transitions that do sometimes exhibit
     * orientation bias filter, sill assumes single context.
     */
    @Test
    public void testMoreComplicatedCase() {
        // use a more lenient epsilon as these cases are harder TODO: lower this as we learn theta and f
        final double epsilon = 0.5;
        final String refContext = "AGT";

        final List<AltSiteRecord> altDeisgnMatrix = new ArrayList<>();
        final Histogram<Integer> refSiteHistogram = new Histogram<>("depth", refContext);
        refSiteHistogram.prefillBins(IntStream.rangeClosed(0, MAX_REF_DEPTH).boxed().toArray(Integer[]::new));

        final int numExamplesPerAllele = 10000; // this many examples, ref and alt, per allele
        final int numAltExamples = 10; // for each alt allele we have this many alt sites
        // TODO: want to have -2 here instead of -1 and still detect the artifact. For that, we must learn theta
        // TODO: sample from a binomial distribution and what not
        final int numGToTF1R2Artifact = ALT_DEPTH - 1; // G -> T transversion only, this many alt sites have F1R2 artifact
        final int numGToAF2R1Artifact = ALT_DEPTH - 1; // G -> A transition only, this many alt sites have F2R1 artifact

        // first create the examples for the G -> T transitions, a fraction of which has read orientation bias
        for (Nucleotide allele : REGULAR_BASES) {
            for (int n = 0; n < numExamplesPerAllele; n++) {
                if (n < numAltExamples && allele != Nucleotide.G) {
                    // Give G -> T transition 95% F1R2 artifact
                    // Give G -> A transition 95% F2R1 artifact
                    final int[] depthCounts = new int[REGULAR_BASES.size()];
                    final int[] f1r2Counts = new int[REGULAR_BASES.size()];

                    depthCounts[allele.ordinal()] = ALT_DEPTH;
                    switch (allele) {
                        case T : f1r2Counts[allele.ordinal()] = numGToTF1R2Artifact; break; // 19 out of 20 alt reads is F1R2
                        case A : f1r2Counts[allele.ordinal()] = ALT_DEPTH - numGToAF2R1Artifact; break; // 19 out of 20 alt reads is F2R1
                        case C : f1r2Counts[allele.ordinal()] = ALT_DEPTH / 2; break; // G -> C should be balanced
                        default : throw new UserException(String.format("We should never reach here but got allele %s", allele));
                    }

                    altDeisgnMatrix.add(new AltSiteRecord(contig, position, refContext, depthCounts,
                            f1r2Counts, DEPTH, allele));
                } else {
                    // create hom ref examples
                    refSiteHistogram.increment(DEPTH);
                }
            }
        }

        // To consider: Can a site be both F1R2 and F2R1? What about the whole reverse complement thing?
        LearnHyperparametersEngine engine = new LearnHyperparametersEngine(refSiteHistogram, altDeisgnMatrix,
                LearnHyperparameters.DEFAULT_CONVERGENCE_THRESHOLD, LearnHyperparameters.DEFAULT_MAX_ITERATIONS, logger);
        Hyperparameters hyperparameters = engine.runEMAlgorithm();

        Assert.assertEquals(engine.effectiveCounts[LearnHyperparametersEngine.ArtifactState.F1R2_T.ordinal()], (double) numAltExamples, epsilon);
        Assert.assertEquals(engine.effectiveCounts[ArtifactState.F2R1_A.ordinal()], (double) numAltExamples, epsilon);
        Assert.assertEquals(engine.effectiveCounts[ArtifactState.SOMATIC_HET.ordinal()], (double) numAltExamples, epsilon);

        // test pi
        Assert.assertEquals(hyperparameters.getPi()[ArtifactState.F1R2_T.ordinal()], (double) numAltExamples/numExamplesPerAllele, 1e-3);
        Assert.assertEquals(hyperparameters.getPi()[ArtifactState.F2R1_A.ordinal()], (double) numAltExamples/numExamplesPerAllele, 1e-3);
        Assert.assertEquals(hyperparameters.getPi()[ArtifactState.SOMATIC_HET.ordinal()], (double) numAltExamples/numExamplesPerAllele, 1e-3);

        // impossible states get 0 probability
        for (ArtifactState z : ArtifactState.getImpossibleStates(Nucleotide.G)){
            Assert.assertEquals(hyperparameters.getPi()[z.ordinal()], 0.0);
        }

    }

    /**
     * Test the case where the alt reads are heavily biased towards F1R2 or F2R1 but not entirely so
     */
    @Test
    public void testEvenMoreComplicated() {

    }

    /**
     * Test multiple contexts, not just the standard old
     */
    @Test
    public void testMultipleContexts() {

    }
}