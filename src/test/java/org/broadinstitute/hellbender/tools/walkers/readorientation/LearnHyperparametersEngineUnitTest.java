package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.ArtifactState;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactType.F1R2;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactType.F2R1;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientationFilterConstants.*;

/**R
 * Created by tsato on 8/14/17.
 */
public class LearnHyperparametersEngineUnitTest extends CommandLineProgramTest {
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
        refSiteHistogram.prefillBins(IntStream.rangeClosed(1, maxDepthForHistograms).boxed().toArray(Integer[]::new));

        for (int n = 0; n < NUM_EXAMPLES; n++){
            // assume 20 in 100 examples is alt at ALT_DEPTH/DEPTH % allele fraction, and alt reads are entirely F1R2
            if (n < NUM_ALT_EXAMPLES){
                altDeisgnMatrix.add(new AltSiteRecord(contig, position, refContext, altSiteDepthCounts,
                        biasedAltF1R2DepthCounts, DEPTH, Nucleotide.T));
            } else {
                refSiteHistogram.increment(DEPTH);
            }
        }

        LearnHyperparametersEngine engine = new LearnHyperparametersEngine(refSiteHistogram, Collections.emptyList(), altDeisgnMatrix,
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
        refSiteHistogram.prefillBins(IntStream.rangeClosed(0, maxDepthForHistograms).boxed().toArray(Integer[]::new));

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
        LearnHyperparametersEngine engine = new LearnHyperparametersEngine(refSiteHistogram, Collections.emptyList(), altDeisgnMatrix,
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

    @Test
    public void testReverseComplement() throws IOException {
        final File altMatrixOutput = File.createTempFile("alt-table", ".tsv");
        final File refHistogramOutput = createTempFile("ref-histogram", "metrics");

        final MetricsFile<?, Integer> refMetricsFile = new MetricsFile<>();

        final AltSiteRecord.AltSiteRecordTableWriter altTableWriter = new AltSiteRecord.AltSiteRecordTableWriter(altMatrixOutput);

        final double epsilon = 1e-3;
        final int numRefExamples = 10000;
        final int refDepth = 150;
        final int numAltExamples = 100;
        final int altDepth = 20;

        final int expectedNumUniqueContexts = 4;
        final List<Triple<String, Nucleotide, ArtifactType>> transitions = Arrays.asList(
                new ImmutableTriple<>("TCT", Nucleotide.G, F1R2), // equivalent to AGA->C F2R1
                new ImmutableTriple<>("AGA", Nucleotide.C, F1R2),
                new ImmutableTriple<>("CAT", Nucleotide.T, F2R1), // equivalent to ATG->A F1R2
                new ImmutableTriple<>("ATG", Nucleotide.A, F1R2),
                new ImmutableTriple<>("TGG", Nucleotide.C, F1R2), // should be recorded as CCA->G F2R1
                new ImmutableTriple<>("GAT", Nucleotide.G, F2R1)); // should be recorded as ATC->C F1R2

        final List<AltSiteRecord> altDesignMatrix = new ArrayList<>();

        for (final Triple<String, Nucleotide, ArtifactType> transition : transitions) {
            final String refContext = transition.getLeft();
            final Nucleotide altAllele = transition.getMiddle();
            final ArtifactType type = transition.getRight();
            final Nucleotide refAllele = Nucleotide.valueOf(refContext.substring(MIDDLE_INDEX, MIDDLE_INDEX + 1));

            final Histogram<Integer> refSiteHistogram = new Histogram<>("depth", refContext);
            refSiteHistogram.prefillBins(bins);

            refSiteHistogram.increment(refDepth, numRefExamples);
            refMetricsFile.addHistogram(refSiteHistogram);

            final int[] depthCounts = new int[REGULAR_BASES.size()];
            final int[] f1r2Counts = new int[REGULAR_BASES.size()];

            depthCounts[altAllele.ordinal()] = altDepth;
            depthCounts[refAllele.ordinal()] = refDepth;

            f1r2Counts[altAllele.ordinal()] = type == F1R2 ? altDepth : 0;
            f1r2Counts[refAllele.ordinal()] = refDepth / 2;

            altTableWriter.writeRecord(new AltSiteRecord(contig, position, refContext, depthCounts,
                    f1r2Counts, altDepth + refDepth, altAllele));
            IntStream.range(0, numAltExamples).forEach(i ->
                    altDesignMatrix.add(new AltSiteRecord(contig, position, refContext, depthCounts,
                            f1r2Counts, altDepth + refDepth, altAllele)));
        }

        refMetricsFile.write(refHistogramOutput);
        altTableWriter.writeAllRecords(altDesignMatrix);
        altTableWriter.close();


        final File hyperparameterOutput = createTempFile("hyperparameters", ".tsv");
        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                        "-alt-table", altMatrixOutput.getAbsolutePath(),
                        "-ref-table", refHistogramOutput.getAbsolutePath(),
                        "-O", hyperparameterOutput.getAbsolutePath()),
                LearnHyperparameters.class.getSimpleName()));

        List<Hyperparameters> priors = Hyperparameters.readHyperparameters(hyperparameterOutput);

        Assert.assertEquals(priors.size(), expectedNumUniqueContexts);
        final Hyperparameters hypForAGA = priors.stream()
                .filter(p -> p.getReferenceContext().equals("AGA"))
                .findFirst()
                .get();

        final Hyperparameters hypForATG = priors.stream()
                .filter(p -> p.getReferenceContext().equals("ATG"))
                .findFirst()
                .get();

        final Hyperparameters hypForCCA = priors.stream()
                .filter(p -> p.getReferenceContext().equals("CCA"))
                .findFirst()
                .get();

        final Hyperparameters hypForATC = priors.stream()
                .filter(p -> p.getReferenceContext().equals("ATC"))
                .findFirst()
                .get();

        if (altTableWriter != null) {
            try {
                altTableWriter.close();
            } catch (IOException e) {
                throw new UserException("Encountered an IO exception while closing the alt table writer", e);
            }
        }

        final double numExamples = numAltExamples + numRefExamples;
        // AGA should get transition to C in both F1R2 and F2R1
        Assert.assertEquals(hypForAGA.getPi(ArtifactState.F1R2_C), numAltExamples / (2 * numExamples), epsilon);
        Assert.assertEquals(hypForAGA.getPi(ArtifactState.F2R1_C), numAltExamples / (2 * numExamples), epsilon);
        Assert.assertEquals(hypForAGA.getPi(ArtifactState.HOM_REF), 2 * numRefExamples / (2 * numExamples), epsilon);

        // ATG should only have one transition to A in F1R2
        Assert.assertEquals(hypForATG.getPi(ArtifactState.F1R2_A), 2 * numAltExamples / (2 * numExamples), epsilon);
        Assert.assertEquals(hypForATG.getPi(ArtifactState.HOM_REF), 2 * numRefExamples / (2 * numExamples), epsilon);

        Assert.assertEquals(hypForCCA.getPi(ArtifactState.F2R1_G), numAltExamples / numExamples, epsilon);
        Assert.assertEquals(hypForATC.getPi(ArtifactState.F1R2_C), numAltExamples / numExamples, epsilon);
    }

    @Test
    public void testMergeHistograms(){
        final int numRefExamples = 10000;
        final int numAltExamples = 100;
        final double epsilon = 1e-3;

        final List<Triple<String, Nucleotide, ArtifactType>> transitions = Arrays.asList(
                new ImmutableTriple<>("TCT", Nucleotide.G, F1R2), // equivalent to AGA->C F2R1
                new ImmutableTriple<>("AGA", Nucleotide.C, F1R2),
                new ImmutableTriple<>("CAT", Nucleotide.T, F2R1), // equivalent to ATG->A F1R2
                new ImmutableTriple<>("ATG", Nucleotide.A, F1R2),
                new ImmutableTriple<>("TGG", Nucleotide.C, F1R2), // should be recorded as CCA->G F2R1
                new ImmutableTriple<>("GAT", Nucleotide.G, F2R1)); // should be recorded as ATC->C F1R2

        // Test Merging Reference Histograms
        final Histogram<Integer> refF1R2TCT = createRefHistograms("TCT");

        final Histogram<Integer> refF1R2AGA = createRefHistograms("AGA");

        final Histogram<Integer> combinedRefAGA = LearnHyperparameters.combineRefHistogramWithRC(
                "AGA", refF1R2AGA, refF1R2TCT);

        Assert.assertEquals((int) combinedRefAGA.get(150).getValue(), 2*numRefExamples);
        Assert.assertEquals((int) combinedRefAGA.getSumOfValues(), 2*numRefExamples);

        // Test Merging Alt Histograms i.e. histograms for alt sites with alt depth = 1
        final List<Histogram<Integer>> altTCT = createAltHistograms("TCT");
        final List<Histogram<Integer>> altAGA = createAltHistograms("AGA");

        final List<Histogram<Integer>> combinedAltAGA = LearnHyperparameters.combineAltHistogramWithRC(altAGA, altTCT);

        // Feeding the input the other way should fail
        try {
            LearnHyperparameters.combineAltHistogramWithRC(altTCT, altAGA);
        } catch (IllegalArgumentException e){
            // Good
        }

        Assert.assertEquals(combinedAltAGA.size(), numSubHistograms);
        Assert.assertEquals(combinedAltAGA.stream()
                .filter(h -> labelToContext(h.getValueLabel()).getLeft().equals("AGA"))
                .count(), numSubHistograms);
        combinedAltAGA.stream().forEach(h ->
                Assert.assertEquals((int) h.getSumOfValues(), 2*100));
        combinedAltAGA.stream().forEach(h ->
                Assert.assertEquals((int) h.get(150).getValue(), 2*100));
    }


    @Test
    public void testMergeDesignMatrices(){
        final int numRefExamples = 10000;
        final int refDepth = 150;
        final int numAltExamples = 100;
        final int altDepth = 20;

        final int expectedNumUniqueContexts = 4;
        final List<Triple<String, Nucleotide, ArtifactType>> transitions1a =
                Arrays.asList(new ImmutableTriple<>("AGA", Nucleotide.C, F1R2));
        final List<Triple<String, Nucleotide, ArtifactType>> transitions1b =
                Arrays.asList(new ImmutableTriple<>("TCT", Nucleotide.G, F1R2));
        final List<Triple<String, Nucleotide, ArtifactType>> transitions1c =
                Arrays.asList(new ImmutableTriple<>("TGG", Nucleotide.C, F1R2));

        final List<AltSiteRecord> altDesignMatrix1a = createDesignMatrix(transitions1a);
        final List<AltSiteRecord> altDesignMatrix1b = createDesignMatrix(transitions1b);
        final List<AltSiteRecord> altDesignMatrix1c = createDesignMatrix(transitions1c);

        // Should be all AGA->C
        LearnHyperparameters.mergeDesignMatrices(altDesignMatrix1a, altDesignMatrix1b);
        Assert.assertEquals(altDesignMatrix1a.stream()
                        .filter(a -> a.getReferenceContext().equals("AGA"))
                        .count(), 2*numAltExamples);
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> a.getAltAllele() == Nucleotide.C)
                .count(), 2*numAltExamples);
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> Arrays.equals(a.getBaseCounts(), new int[]{0,20,150,0}))
                .count(), 2*numAltExamples);
        // Half of the sites should be heavily biased for F1R2 with the other half biased for F2R1
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> Arrays.equals(a.getF1R2Counts(), new int[]{0,0,75,0}))
                .count(), numAltExamples);
        Assert.assertEquals(altDesignMatrix1a.stream()
                .filter(a -> Arrays.equals(a.getF1R2Counts(), new int[]{0,20,75,0}))
                .count(), numAltExamples);


        try {
            // Merging the wrong direction should throw an error
            LearnHyperparameters.mergeDesignMatrices(altDesignMatrix1b, altDesignMatrix1a);
            Assert.fail();
        } catch (IllegalArgumentException e){
            // Good
        }

        try {
            // Merging non-matching contexts shoudl throw and error
            LearnHyperparameters.mergeDesignMatrices(altDesignMatrix1a, altDesignMatrix1c);
            Assert.fail();
        } catch (IllegalArgumentException e){
            // Good
        }

    }

    private List<AltSiteRecord> createDesignMatrix(final List<Triple<String, Nucleotide, ArtifactType>> transitions) {
        final int altDepth = 20;
        final int refDepth = 150;
        final int numAltExamples = 100;

        final List<AltSiteRecord> altDesignMatrix = new ArrayList<>(transitions.size()*numAltExamples);
        for (final Triple<String, Nucleotide, ArtifactType> transition : transitions) {
            final String refContext = transition.getLeft();
            final Nucleotide altAllele = transition.getMiddle();
            final ArtifactType type = transition.getRight();
            final Nucleotide refAllele = Nucleotide.valueOf(refContext.substring(MIDDLE_INDEX, MIDDLE_INDEX + 1));

            final int[] depthCounts = new int[REGULAR_BASES.size()];
            final int[] f1r2Counts = new int[REGULAR_BASES.size()];

            depthCounts[altAllele.ordinal()] = altDepth;
            depthCounts[refAllele.ordinal()] = refDepth;

            f1r2Counts[altAllele.ordinal()] = type == F1R2 ? altDepth : 0;
            f1r2Counts[refAllele.ordinal()] = refDepth / 2;

            // altTableWriter.writeRecord(new AltSiteRecord(contig, position, refContext, depthCounts,
            // f1r2Counts, altDepth + refDepth, altAllele));
            IntStream.range(0, numAltExamples).forEach(i ->
                    altDesignMatrix.add(new AltSiteRecord(contig, position, refContext, depthCounts,
                            f1r2Counts, altDepth + refDepth, altAllele)));
        }

        return altDesignMatrix;
    }

    // Returns <ref, alt>
    private Histogram<Integer> createRefHistograms(String refContext) {
        int refDepth = 150;
        int numRefExamples = 10000;

        final Histogram<Integer> refSiteHistogram = new Histogram<>(binName, refContext);
        refSiteHistogram.prefillBins(bins);

        refSiteHistogram.increment(refDepth, numRefExamples);

        return refSiteHistogram;
    }

    private List<Histogram<Integer>> createAltHistograms(final String refContext) {
        int refDepth = 150;
        int numAltExamples = 100;

        final List<Histogram<Integer>> altComputationalHistograms = new ArrayList<>(
                ArtifactType.values().length * (REGULAR_BASES.size() - 1));
        for (Nucleotide altAllele : REGULAR_BASES){
            if (altAllele.toString().equals(refContext.substring(MIDDLE_INDEX, MIDDLE_INDEX+1))){
                continue;
            }

            for (ArtifactType type : ArtifactType.values()){
                final String columnLabel = contextToLabel(refContext, altAllele, type);
                final Histogram<Integer> altComputationalHistogram = new Histogram<>(binName, columnLabel);
                altComputationalHistogram.prefillBins(bins);
                altComputationalHistogram.increment(refDepth, numAltExamples);
                altComputationalHistograms.add(altComputationalHistogram);
            }
        }

        return altComputationalHistograms;
    }
}