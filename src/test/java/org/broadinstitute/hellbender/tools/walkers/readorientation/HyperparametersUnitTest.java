package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.NUM_STATES;

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
}