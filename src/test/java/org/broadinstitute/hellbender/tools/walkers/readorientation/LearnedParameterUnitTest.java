package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Optional;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2FilterConstants.NUM_STATES;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactState.*;

public class LearnedParameterUnitTest {
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

        final String referenceContext1 = "ACT";
        final String referenceContext2 = "GTA";
        final String referenceContext3 = "CCC";

        final LearnedParameterCollection learnedParameterCollectionBefore = new LearnedParameterCollection("sample");
        learnedParameterCollectionBefore.set(new LearnedParameter(referenceContext1, pi, numExamples1, numAltExamples1));
        learnedParameterCollectionBefore.set(new LearnedParameter(referenceContext2, pi, numExamples2, numAltExamples2));
        learnedParameterCollectionBefore.set(new LearnedParameter(referenceContext3, pi, numExamples3, numAltExamples3));

        final File table = File.createTempFile("prior", ".tsv");
        learnedParameterCollectionBefore.writeArtifactPriors(table, ParameterType.ARTIFACT_PRIOR);
        final LearnedParameterCollection learnedParameterCollectionAfter = LearnedParameterCollection.readArtifactPriors(table, ParameterType.ARTIFACT_PRIOR);

        final LearnedParameter ap1 = learnedParameterCollectionAfter.get(referenceContext1).get();
        ArrayAsserts.assertArrayEquals(ap1.getParameters(), pi, EPSILON);
        Assert.assertEquals(ap1.getNumExamples(), numExamples1);
        Assert.assertEquals(ap1.getNumAltExamples(), numAltExamples1);

        LearnedParameter ap2 = learnedParameterCollectionAfter.get(referenceContext2).get();
        ArrayAsserts.assertArrayEquals(ap2.getParameters(), pi, EPSILON);
        Assert.assertEquals(ap2.getNumExamples(), numExamples2);
        Assert.assertEquals(ap2.getNumAltExamples(), numAltExamples2);

        LearnedParameter ap3 = learnedParameterCollectionAfter.get(referenceContext3).get();
        ArrayAsserts.assertArrayEquals(ap3.getParameters(), pi, EPSILON);
        Assert.assertEquals(ap3.getNumExamples(), numExamples3);
        Assert.assertEquals(ap3.getNumAltExamples(), numAltExamples3);
    }

    @Test
    public void testRevComp() throws IOException {
        final int numExamples = 1000;
        final int numAltExamples = 10;
        final String referenceContext1 = "ACT";
        final String referenceContext2 = "CCC";
        final double[] pi1 = new double[F1R2FilterConstants.NUM_STATES];
        ArtifactState.setStatePrior(pi1, 0.3, F1R2_A);
        ArtifactState.setStatePrior(pi1, 0.7, HOM_REF);

        final double[] pi2 = new double[F1R2FilterConstants.NUM_STATES];
        ArtifactState.setStatePrior(pi2, 0.05, F2R1_C);
        ArtifactState.setStatePrior(pi2, 0.8, HOM_REF);
        ArtifactState.setStatePrior(pi2, 0.15, SOMATIC_HET);

        final LearnedParameterCollection learnedParameterCollectionBefore = new LearnedParameterCollection("sample");
        learnedParameterCollectionBefore.set(new LearnedParameter(referenceContext1, pi1, numExamples, numAltExamples));
        learnedParameterCollectionBefore.set(new LearnedParameter(referenceContext2, pi2, numExamples, numAltExamples));

        final File table = File.createTempFile("prior", ".tsv");
        learnedParameterCollectionBefore.writeArtifactPriors(table, ParameterType.ARTIFACT_PRIOR);

        final LearnedParameterCollection learnedParameterCollectionAfter = LearnedParameterCollection.readArtifactPriors(table, ParameterType.ARTIFACT_PRIOR);
        final LearnedParameter ap1 = learnedParameterCollectionAfter.get(referenceContext1).get();
        final LearnedParameter ap1rc = learnedParameterCollectionAfter.get(SequenceUtil.reverseComplement(referenceContext1)).get();

        Assert.assertEquals(ap1.getReferenceContext(), SequenceUtil.reverseComplement(ap1rc.getReferenceContext()));
        for (ArtifactState state : ArtifactState.values()){
            Assert.assertEquals(ap1.getParameter(state), ap1rc.getParameter(state.getRevCompState()));
        }

        final LearnedParameter ap2 = learnedParameterCollectionAfter.get(referenceContext2).get();
        final LearnedParameter ap2rc = learnedParameterCollectionAfter.get(SequenceUtil.reverseComplement(referenceContext2)).get();

        Assert.assertEquals(ap2.getReferenceContext(), SequenceUtil.reverseComplement(ap2rc.getReferenceContext()));
        for (ArtifactState state : ArtifactState.values()){
            Assert.assertEquals(ap2.getParameter(state), ap2rc.getParameter(state.getRevCompState()));
        }

    }

    @DataProvider(name = "dataForGet")
    public Object[][] dataForGet() {
        final double[] pi = new double[F1R2FilterConstants.NUM_STATES];
        ArtifactState.setStatePrior(pi, 0.2, F1R2_A);
        ArtifactState.setStatePrior(pi, 0.8, HOM_REF);

        final double[] expectedPiRC = new double[F1R2FilterConstants.NUM_STATES];
        ArtifactState.setStatePrior(expectedPiRC, 0.8, HOM_REF);
        ArtifactState.setStatePrior(expectedPiRC, 0.2, F2R1_T);

        final double[] pi2 = new double[F1R2FilterConstants.NUM_STATES];
        ArtifactState.setStatePrior(pi2, 0.1, F1R2_A);
        ArtifactState.setStatePrior(pi2, 0.2, F1R2_G);
        ArtifactState.setStatePrior(pi2, 0.3, F2R1_C);
        ArtifactState.setStatePrior(pi2, 0.4, HOM_REF);

        final double[] expectedPiRC2 = new double[F1R2FilterConstants.NUM_STATES];
        ArtifactState.setStatePrior(expectedPiRC2, 0.1, F2R1_T);
        ArtifactState.setStatePrior(expectedPiRC2, 0.2, F2R1_C);
        ArtifactState.setStatePrior(expectedPiRC2, 0.3, F1R2_G);
        ArtifactState.setStatePrior(expectedPiRC2, 0.4, HOM_REF);

        return new Object[][]{
                {pi, expectedPiRC},
                {pi2, expectedPiRC2 }};
    }

    @Test(dataProvider = "dataForGet")
    public void testGet(final double[] pi, final double[] expectedRevCompPi) {
        final double epsilon = 1e-6;

        final String context = "ATG";
        final String revComp = "CAT";

        final LearnedParameterCollection priors = new LearnedParameterCollection("sample");
        priors.set(new LearnedParameter(context, pi, 100, 10));
        Optional<LearnedParameter> priorRevComp = priors.get(revComp);

        ArrayAsserts.assertArrayEquals(expectedRevCompPi, priorRevComp.get().getParameters(), epsilon);
    }
}