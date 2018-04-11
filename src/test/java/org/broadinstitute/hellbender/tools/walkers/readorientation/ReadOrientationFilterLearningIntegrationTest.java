package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.walkers.annotator.ReferenceBases.REFERENCE_BASES_KEY;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactType.F1R2;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactType.F2R1;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.ArtifactState.*;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.ArtifactState;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientationFilterConstants.CANONICAL_KMERS;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.*;

/**
 * Created by tsato on 8/1/17.
 */
public class ReadOrientationFilterLearningIntegrationTest extends CommandLineProgramTest {
    /**
     * Test the tool on a real bam to make sure that it does not crash
     */
    @Test
    public void testOnRealBam() throws IOException {
        final File refMetrics = createTempFile("ref", ".table");
        final File altMetrics = createTempFile("alt", ".table");
        final File altTable = createTempFile("alt", ".table");

        final String testDir = toolsTestDir + "read_orientation_filter/";
        // final String hapmapBamSnippet = testDir + "hapmap-20-plex-chr20-ROF.bam";
        final String hapmapBamSnippet = testDir + "hapmap-20-plex-chr-20-21-read-orientation.bam";

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                        "-R", b37_reference_20_21,
                        "-I", hapmapBamSnippet,
                        "-" + CollectDataForReadOrientationFilter.ALT_DATA_TABLE_SHORT_NAME, altTable.getAbsolutePath(),
                        "-" + CollectDataForReadOrientationFilter.REF_SITE_METRICS_SHORT_NAME, refMetrics.getAbsolutePath(),
                        "-alt-hist", altMetrics.getAbsolutePath()),
                CollectDataForReadOrientationFilter.class.getSimpleName()));

        int lineCount = (int) Files.lines(Paths.get(refMetrics.getAbsolutePath())).filter(l -> l.matches("^[0-9].+")).count();

        // Ensure that we print every bin, even when the count is 0
        Assert.assertEquals(lineCount, ReadOrientationFilterConstants.maxDepthForHistograms);

        // Run LearnHyperparameter
        final File hyperparameters = createTempFile("hyperparameters", ".tsv");
        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                    "-alt-table", altTable.getAbsolutePath(),
                    "-alt-histogram",  altMetrics.getAbsolutePath(),
                    "-ref-table", refMetrics.getAbsolutePath(),
                    "-O", hyperparameters.getAbsolutePath()),
                LearnHyperparameters.class.getSimpleName()));

        final List<Hyperparameters> hyperparametersList = Hyperparameters.readHyperparameters(hyperparameters);
        final Hyperparameters hypsForACT = Hyperparameters.searchByContext(hyperparametersList, "ACT").get();
        final Hyperparameters hypsForAGT = Hyperparameters.searchByContext(hyperparametersList, "AGT").get();

        // NOT TRUE! FIX
        ArrayAsserts.assertArrayEquals(hypsForACT.getPi(), hypsForAGT.getPi(), 1e-5);

        // Run Mutect 2
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                    "-I", hapmapBamSnippet,
                    "-" + M2ArgumentCollection.TUMOR_SAMPLE_SHORT_NAME, "SM-CEMAH",
                    "-R", b37_reference_20_21,
                    "--table", hyperparameters.getAbsolutePath(),
                    "-O", unfilteredVcf.getAbsolutePath()),
                Mutect2.class.getSimpleName()));

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                        "-V", unfilteredVcf.getAbsolutePath(),
                        "-R", b37_reference_20_21,
                        "-O", filteredVcf.getAbsolutePath()),
                FilterMutectCalls.class.getSimpleName()));


        // These artifacts have been verified manually
        // The triple is of type (Position, Artifact Type, Expected Source of Prior Probability)
        // Prior for e.g. TGA->A F2R1 should come from TCA->T F1R2
        final List<Triple<Integer, ArtifactType, ArtifactState>> knownArtifacts = Arrays.asList(
                new ImmutableTriple<>(5296933, F1R2, F1R2_T), // CGT->T, which is not canonical
                new ImmutableTriple<>(23421079, F2R1, F2R1_G), // CAC->G
                new ImmutableTriple<>(34144749, F2R1, F2R1_A), // TGA->A, which is not canonical. Expect F1R2_T
                new ImmutableTriple<>(62165528, F1R2, F1R2_A)); // CGC->A

        for (final Triple<Integer, ArtifactType, ArtifactState> artifact : knownArtifacts) {
            final int position = artifact.getLeft();
            final ArtifactType type = artifact.getMiddle();
            final ArtifactState expectedSourceOfPrior = artifact.getRight();

            Optional<VariantContext> variant = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                    .filter(vc -> vc.getStart() == position).findFirst();
            Assert.assertTrue(variant.isPresent());

            // Check that the correct prior was added to the format field by Mutect
            final double prior = GATKProtectedVariantContextUtils.getAttributeAsDouble(variant.get().getGenotype(0), ROF_PRIOR_KEY, -1.0);
            final String refContext = variant.get().getAttributeAsString(REFERENCE_BASES_KEY, "");
            final Hyperparameters hyps = Hyperparameters.searchByContext(hyperparametersList, refContext).get();

            Assert.assertEquals(prior, hyps.getPi(expectedSourceOfPrior), 1e-3);

            final String expectedFilter = type == F1R2 ? F1R2_ARTIFACT_FILTER_NAME : F2R1_ARTIFACT_FILTER_NAME;

            // Check that the expected filters were applied
            Assert.assertTrue(variant.get().getFilters().contains(expectedFilter));
            Assert.assertEquals(GATKProtectedVariantContextUtils.getAttributeAsString(variant.get().getGenotype(0), ROF_TYPE_KEY, null),
                    type.toString());
        }
    }
}