package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.NUM_STATES;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.ArtifactState;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.ArtifactState.*;

/**
 * Created by tsato on 9/6/17.
 */
public class HyperparametersUnitTest {
    private static final double EPSILON = 1e-4;

    @Test
    public void testIO() throws IOException {
        final double[] pi = new double[NUM_STATES];
        Arrays.fill(pi, 1.0/NUM_STATES);

        final int numExamples1 = 30;
        final int numExamples2 = 60;
        final int numExamples3 = 90;

        final int numAltExamples1 = 10;
        final int numAltExamples2 = 50;
        final int numAltExamples3 = 1;



        List<Hyperparameters> hyperparametersBefore = Arrays.asList(
                new Hyperparameters("ACT", pi, numExamples1, numAltExamples1),
                new Hyperparameters("GTA", pi, numExamples2, numAltExamples2),
                new Hyperparameters("CCC", pi, numExamples3, numAltExamples3));
        final File table = File.createTempFile("hyperparameters", ".tsv");

        Hyperparameters.writeHyperparameters(hyperparametersBefore, table);
        List<Hyperparameters> hyperparametersAfter = Hyperparameters.readHyperparameters(table);

        Hyperparameters hyp1 = hyperparametersAfter.get(0);
        ArrayAsserts.assertArrayEquals(hyp1.getPi(), pi, EPSILON);
        Assert.assertEquals(hyp1.getNumExamples(), numExamples1);
        Assert.assertEquals(hyp1.getNumAltExamples(), numAltExamples1);

        Hyperparameters hyp2 = hyperparametersAfter.get(1);
        ArrayAsserts.assertArrayEquals(hyp2.getPi(), pi, EPSILON);
        Assert.assertEquals(hyp2.getNumExamples(), numExamples2);
        Assert.assertEquals(hyp2.getNumAltExamples(), numAltExamples2);

        Hyperparameters hyp3 = hyperparametersAfter.get(2);
        ArrayAsserts.assertArrayEquals(hyp3.getPi(), pi, EPSILON);
        Assert.assertEquals(hyp3.getNumExamples(), numExamples3);
        Assert.assertEquals(hyp3.getNumAltExamples(), numAltExamples3);
    }

    @Test
    public void testRevComp() throws IOException {
        final int numExamples = 1000;
        final int numAltExamples = 10;
        final double[] pi = new double[LearnHyperparametersEngine.ArtifactState.values().length];
        final String context1 = "ACT";
        final String context2 = "CCC";
        List<Hyperparameters> hyperparametersBefore = Arrays.asList(
                new Hyperparameters(context1, pi, numExamples, numAltExamples),
                new Hyperparameters(context2, pi, numExamples, numAltExamples));
        final File table = File.createTempFile("hyperparameters", ".tsv");

        Hyperparameters.writeHyperparameters(hyperparametersBefore, table);
        List<Hyperparameters> hyperparametersAfter = Hyperparameters.readHyperparameters(table);
        Hyperparameters hypForContext1 = Hyperparameters.searchByContext(hyperparametersAfter, context1).get();
        Hyperparameters hypForContext1RevComp = Hyperparameters.searchByContext(hyperparametersAfter,
                SequenceUtil.reverseComplement(context1)).get();
        Assert.assertEquals(hypForContext1.getReferenceContext(), hypForContext1RevComp.getReferenceContext());

        Hyperparameters hypForContext2 = Hyperparameters.searchByContext(hyperparametersAfter, context2).get();
        Hyperparameters hypForContext2RevComp = Hyperparameters.searchByContext(hyperparametersAfter,
                SequenceUtil.reverseComplement(context2)).get();
        Assert.assertEquals(hypForContext2.getReferenceContext(), hypForContext2RevComp.getReferenceContext());
    }

    @DataProvider(name = "searchByContextData")
    public Object[][] searchByContextData() {
        final double[] pi = new double[ArtifactState.values().length];
        pi[F1R2_A.ordinal()] = 0.2;
        pi[HOM_REF.ordinal()] = 0.8;

        final double[] expectedPiRC = new double[ArtifactState.values().length];
        expectedPiRC[HOM_REF.ordinal()] = 0.8;
        expectedPiRC[F2R1_T.ordinal()] = 0.2;

        final double[] pi2 = new double[ArtifactState.values().length];
        pi2[F1R2_A.ordinal()] = 0.1;
        pi2[F1R2_G.ordinal()] = 0.2;
        pi2[F2R1_C.ordinal()] = 0.3;
        pi2[HOM_REF.ordinal()] = 0.4;

        final double[] expectedPiRC2 = new double[ArtifactState.values().length];
        expectedPiRC2[F2R1_T.ordinal()] = 0.1;
        expectedPiRC2[F2R1_C.ordinal()] = 0.2;
        expectedPiRC2[F1R2_G.ordinal()] = 0.3;
        expectedPiRC2[HOM_REF.ordinal()] = 0.4;

        return new Object[][]{
                {pi, expectedPiRC},
                {pi2, expectedPiRC2 }};
    }

    @Test(dataProvider = "searchByContextData")
    public void testSearchByContext(final double[] pi, final double[] expectedRevCompPi) {
        final double epsilon = 1e-6;

        final String context = "ATG";
        final String revComp = SequenceUtil.reverseComplement(context);

        final Hyperparameters hyps = new Hyperparameters(context, pi, 100, 10);
        Optional<Hyperparameters> h = Hyperparameters.searchByContext(Arrays.asList(hyps), context);
        Optional<Hyperparameters> hRC = Hyperparameters.searchByContext(Arrays.asList(hyps), revComp);

        ArrayAsserts.assertArrayEquals(expectedRevCompPi, hRC.get().getPi(), epsilon);
    }

}