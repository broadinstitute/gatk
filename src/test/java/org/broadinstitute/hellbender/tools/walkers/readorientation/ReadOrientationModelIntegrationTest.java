package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReferenceBases;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.StreamSupport;

public class ReadOrientationModelIntegrationTest extends CommandLineProgramTest {
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
                        "--" + CollectF1R2Counts.ALT_DATA_TABLE_LONG_NAME, altTable.getAbsolutePath(),
                        "--" + CollectF1R2Counts.REF_SITE_METRICS_LONG_NAME, refMetrics.getAbsolutePath(),
                        "--" + CollectF1R2Counts.ALT_DEPTH1_HISTOGRAM_LONG_NAME, altMetrics.getAbsolutePath()),
                CollectF1R2Counts.class.getSimpleName()));

        int lineCount = (int) Files.lines(Paths.get(refMetrics.getAbsolutePath())).filter(l -> l.matches("^[0-9].+")).count();

        // Ensure that we print every bin, even when the count is 0
        Assert.assertEquals(lineCount, F1R2FilterConstants.DEFAULT_MAX_DEPTH);

        // Run the prior probabilities of artifact
        final File priorTable = createTempFile("prior", ".tsv");
        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                    "--" + CollectF1R2Counts.ALT_DATA_TABLE_LONG_NAME, altTable.getAbsolutePath(),
                    "--" + CollectF1R2Counts.ALT_DEPTH1_HISTOGRAM_LONG_NAME,  altMetrics.getAbsolutePath(),
                    "--" + CollectF1R2Counts.REF_SITE_METRICS_LONG_NAME, refMetrics.getAbsolutePath(),
                    "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, priorTable.getAbsolutePath()),
                LearnReadOrientationModel.class.getSimpleName()));

        final ArtifactPriorCollection artifactPriorCollection = ArtifactPriorCollection.readArtifactPriors(priorTable);

        // Run Mutect 2
        final File unfilteredVcf = GATKBaseTest.createTempFile("unfiltered", ".vcf");
        final File filteredVcf = GATKBaseTest.createTempFile("filtered", ".vcf");
        final File bamout = GATKBaseTest.createTempFile("SM-CEMAH", ".bam");

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                    "-I", hapmapBamSnippet,
                    "-R", b37_reference_20_21,
                    "--" + M2ArgumentCollection.ARTIFACT_PRIOR_TABLE_NAME, priorTable.getAbsolutePath(),
                    "-O", unfilteredVcf.getAbsolutePath(),
                    "-bamout", bamout.getAbsolutePath()),
                Mutect2.class.getSimpleName()));

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                        "-V", unfilteredVcf.getAbsolutePath(),
                        "-R", b37_reference_20_21,
                        "-O", filteredVcf.getAbsolutePath(),
                        "--" + M2FiltersArgumentCollection.FALSE_DISCOVERY_RATE_LONG_NAME, "0.04",
                FilterMutectCalls.class.getSimpleName())));


        // These artifacts have been verified manually
        // The pair is of type (Position, Expected Source of Prior Probability)
        // Prior for e.g. TGA->A F2R1 should come from TCA->T F1R2
        final List<Triple<Integer, ReadOrientation, ArtifactState>> knownArtifacts = Arrays.asList(
                new ImmutableTriple<>(23421079, ReadOrientation.F2R1, ArtifactState.F2R1_G), // CAC->G F2R1
                new ImmutableTriple<>(34144749, ReadOrientation.F2R1, ArtifactState.F2R1_A), // TGA->A F2R1, ref context is not canonical
                new ImmutableTriple<>(62165528, ReadOrientation.F1R2, ArtifactState.F1R2_A)); // CGC->A F1R2

        for (final Triple<Integer, ReadOrientation, ArtifactState> artifact : knownArtifacts) {
            final int position = artifact.getLeft();
            final ReadOrientation expectedReadOrientaiton = artifact.getMiddle();
            final ArtifactState expectedSourceOfPrior = artifact.getRight();

            Optional<VariantContext> variant = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                    .filter(vc -> vc.getStart() == position).findFirst();
            Assert.assertTrue(variant.isPresent());

            // Check that the correct prior was added to the format field by Mutect
            final double prior = GATKProtectedVariantContextUtils.getAttributeAsDouble(variant.get().getGenotype(0), GATKVCFConstants.ROF_PRIOR_KEY, -1.0);
            final String refBases = variant.get().getAttributeAsString(ReferenceBases.REFERENCE_BASES_KEY, "");
            final String refContext = ReferenceBases.getNMiddleBases(refBases, F1R2FilterConstants.REFERENCE_CONTEXT_SIZE);
            final ArtifactPrior ap = artifactPriorCollection.get(refContext).get();

            Assert.assertEquals(prior, ap.getPi(expectedSourceOfPrior), 1e-3);

            // Check that the expected filters were applied
            Assert.assertTrue(variant.get().getFilters().contains(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME));
            Assert.assertEquals(GATKProtectedVariantContextUtils.getAttributeAsString(variant.get().getGenotype(0),
                    GATKVCFConstants.ROF_TYPE_KEY, null), expectedReadOrientaiton.toString());
        }
    }
}