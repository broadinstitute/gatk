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

        final double[] f1 = new double[NUM_STATES];
        Arrays.fill(f1, 0.1);

        final double[] f2 = new double[NUM_STATES];
        Arrays.fill(f2, 0.2);

        final double[] f3 = new double[NUM_STATES];
        Arrays.fill(f3, 0.3);

        final double[] theta1 = new double[NUM_STATES];
        Arrays.fill(theta1, 0.9);

        final double[] theta2 = new double[NUM_STATES];
        Arrays.fill(theta2, 0.8);

        final double[] theta3 = new double[NUM_STATES];
        Arrays.fill(theta3, 0.7);

        final int numExamples1 = 30;
        final int numExamples2 = 60;
        final int numExamples3 = 90;

        final int numAltExamples1 = 10;
        final int numAltExamples2 = 50;
        final int numAltExamples3 = 1;



        List<Hyperparameters> hyperparametersBefore = Arrays.asList(
                new Hyperparameters("ACT", pi, f1, theta1, numExamples1, numAltExamples1),
                new Hyperparameters("GTA", pi, f2, theta2, numExamples2, numAltExamples2),
                new Hyperparameters("CCC", pi, f3, theta3, numExamples3, numAltExamples3));
        final File table = File.createTempFile("hyperparameters", ".tsv");

        Hyperparameters.writeHyperparameters(hyperparametersBefore, table);
        List<Hyperparameters> hyperparametersAfter = Hyperparameters.readHyperparameters(table);

        Hyperparameters hyp1 = hyperparametersAfter.get(0);
        ArrayAsserts.assertArrayEquals(hyp1.getPi(), pi, EPSILON);
        ArrayAsserts.assertArrayEquals(hyp1.getF(), f1, EPSILON);
        ArrayAsserts.assertArrayEquals(hyp1.getTheta(), theta1, EPSILON);
        Assert.assertEquals(hyp1.getNumExamples(), numExamples1);
        Assert.assertEquals(hyp1.getNumAltExamples(), numAltExamples1);

        Hyperparameters hyp2 = hyperparametersAfter.get(1);
        ArrayAsserts.assertArrayEquals(hyp2.getPi(), pi, EPSILON);
        ArrayAsserts.assertArrayEquals(hyp2.getF(), f2, EPSILON);
        ArrayAsserts.assertArrayEquals(hyp2.getTheta(), theta2, EPSILON);
        Assert.assertEquals(hyp2.getNumExamples(), numExamples2);
        Assert.assertEquals(hyp2.getNumAltExamples(), numAltExamples2);

        Hyperparameters hyp3 = hyperparametersAfter.get(2);
        ArrayAsserts.assertArrayEquals(hyp3.getPi(), pi, EPSILON);
        ArrayAsserts.assertArrayEquals(hyp3.getF(), f3, EPSILON);
        ArrayAsserts.assertArrayEquals(hyp3.getTheta(), theta3, EPSILON);
        Assert.assertEquals(hyp3.getNumExamples(), numExamples3);
        Assert.assertEquals(hyp3.getNumAltExamples(), numAltExamples3);
    }
}