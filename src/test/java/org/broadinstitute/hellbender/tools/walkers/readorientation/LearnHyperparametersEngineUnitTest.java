package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.State;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by tsato on 8/14/17.
 */
public class LearnHyperparametersEngineUnitTest {
    private static final double EPSILON = 1e-4;

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
        final RefSiteHistogram refSiteHistogram = new RefSiteHistogram(refContext);

        for (int n = 0; n < NUM_EXAMPLES; n++){
            // assume 20 in 100 examples is alt at ALT_DEPTH/DEPTH % allele fraction, and alt reads are entirely F1R2
            if (n < NUM_ALT_EXAMPLES){
                altDeisgnMatrix.add(new AltSiteRecord(refContext, altSiteDepthCounts, biasedAltF1R2DepthCounts, DEPTH, Nucleotide.T));
            } else {
                refSiteHistogram.increment(DEPTH);
            }
        }

        LearnHyperparametersEngine engine = new LearnHyperparametersEngine(refSiteHistogram, altDeisgnMatrix,
                LearnHyperparameters.DEFAULT_CONVERGENCE_THRESHOLD, LearnHyperparameters.DEFAULT_MAX_ITERATIONS);
        Hyperparameters hyperparameters = engine.runEMAlgorithm(logger);

        Assert.assertEquals(engine.effectiveCounts[State.F1R2_T.ordinal()], (double) NUM_F1R2_EXAMPLES, EPSILON);
        Assert.assertEquals(engine.effectiveCounts[State.HOM_REF.ordinal()], (double) NUM_EXAMPLES - NUM_ALT_EXAMPLES, EPSILON);

        // check pi
        Assert.assertEquals(hyperparameters.getPi()[State.F1R2_T.ordinal()], (double) NUM_ALT_EXAMPLES/NUM_EXAMPLES, EPSILON);
        Assert.assertEquals(hyperparameters.getPi()[State.HOM_REF.ordinal()], 1.0 - (double) NUM_ALT_EXAMPLES/NUM_EXAMPLES, EPSILON);

        // We expect the model to learn the correct allele fraction for the case of an artifact
        Assert.assertEquals(hyperparameters.getF()[State.F1R2_T.ordinal()], (double) ALT_DEPTH /DEPTH, EPSILON);
        Assert.assertEquals(hyperparameters.getF()[State.F2R1_T.ordinal()], 0.0, EPSILON);
        Assert.assertEquals(hyperparameters.getF()[State.GERMLINE_HET.ordinal()], 0.0, EPSILON);
        Assert.assertEquals(hyperparameters.getF()[State.SOMATIC_HET.ordinal()], 0.0, EPSILON);
        Assert.assertEquals(hyperparameters.getF()[State.HOM_REF.ordinal()], 0.0, EPSILON);

        // And theta?
    }


    /**
     * Now test the case where not all of the transitions have orientation bias. And on the transitions that do,
     * not all of the sites have orientation bias. Still assumes single context.
     */
    @Test
    public void testMoreComplicatedCase() {
        // use a more lenient epsilon as these cases are harder TODO: lower as we learn theta, f
        final double epsilon = 0.5;
        final String refContext = "AGT";

        final List<AltSiteRecord> altDeisgnMatrix = new ArrayList<>();
        final RefSiteHistogram refSiteHistogram = new RefSiteHistogram(refContext);

        final int numExamplesPerAllele = 10000; // this many examples, ref and alt, per allele
        final int numAltExamples = 10; // for each alt allele we have this many alt sites
        // TODO: want to have - 2 here instead of -1 and still detect the artifact. For that, we must learn theta
        // TODO: make this test more realistic with sampling from a binomial distribution and what not
        final int numGToTF1R2Artifact = ALT_DEPTH - 1; // G -> T transversion only, this many alt sites have F1R2 artifact
        final int numGToAF2R1Artifact = ALT_DEPTH - 1; // G -> A transition only, this many alt sites have F2R1 artifact

        // first create the examples for the G -> T transitions, a fraction of which has read orientation bias
        List<Nucleotide> alleles = Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T);
        for (Nucleotide allele : alleles) {
            for (int n = 0; n < numExamplesPerAllele; n++) {
                if (n < numAltExamples && allele != Nucleotide.G) {
                    // Give G -> T transition 90% F1R2 artifact
                    // Give G -> A transition 95% F2R1 artifact
                    final int[] depthCounts = new int[alleles.size()];
                    final int[] f1r2Counts = new int[alleles.size()];

                    depthCounts[allele.ordinal()] = ALT_DEPTH;
                    switch (allele) {
                        case T : f1r2Counts[allele.ordinal()] = numGToTF1R2Artifact; break; // 19 out of 20 alt reads is F1R2
                        case A : f1r2Counts[allele.ordinal()] = ALT_DEPTH - numGToAF2R1Artifact; break; // 19 out of 20 alt reads is F2R1
                        case C : f1r2Counts[allele.ordinal()] = ALT_DEPTH / 2; break; // G -> C should be balanced
                        default : throw new UserException(String.format("We should never reach here but got allele %s", allele));
                    }

                    altDeisgnMatrix.add(new AltSiteRecord(refContext, depthCounts, f1r2Counts, DEPTH, allele));
                } else {
                    // create hom ref examples
                    refSiteHistogram.increment(DEPTH);
                }
            }
        }

        // To consider: Can a site be both F1R2 and F2R1? What about the whole reverse complement thing?
        LearnHyperparametersEngine engine = new LearnHyperparametersEngine(refSiteHistogram, altDeisgnMatrix,
                LearnHyperparameters.DEFAULT_CONVERGENCE_THRESHOLD, LearnHyperparameters.DEFAULT_MAX_ITERATIONS);
        Hyperparameters hyperparameters = engine.runEMAlgorithm(logger);

        Assert.assertEquals(engine.effectiveCounts[State.F1R2_T.ordinal()], (double) numAltExamples, epsilon);
        Assert.assertEquals(engine.effectiveCounts[State.F2R1_A.ordinal()], (double) numAltExamples, epsilon);
        Assert.assertEquals(engine.effectiveCounts[State.SOMATIC_HET.ordinal()], (double) numAltExamples, epsilon);

        // test pi
        Assert.assertEquals(hyperparameters.getPi()[State.F1R2_T.ordinal()], (double) numAltExamples/numExamplesPerAllele, 1e-3);
        Assert.assertEquals(hyperparameters.getPi()[State.F2R1_A.ordinal()], (double) numAltExamples/numExamplesPerAllele, 1e-3);
        Assert.assertEquals(hyperparameters.getPi()[State.SOMATIC_HET.ordinal()], (double) numAltExamples/numExamplesPerAllele, 1e-3);

        // impossible states get 0 probability
        for (State z : State.getImpossibleStates(Nucleotide.G)){
            Assert.assertEquals(hyperparameters.getPi()[z.ordinal()], 0.0);
        }

        // test f
        Assert.assertEquals(hyperparameters.getF()[State.F1R2_T.ordinal()], (double) ALT_DEPTH/DEPTH, 1e-3);
        Assert.assertEquals(hyperparameters.getF()[State.F2R1_A.ordinal()], (double) ALT_DEPTH/DEPTH, 1e-3);
        Assert.assertEquals(hyperparameters.getF()[State.SOMATIC_HET.ordinal()], (double) ALT_DEPTH/DEPTH, 1e-3);

        // test theta
        // TODO: implement learning of theta

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