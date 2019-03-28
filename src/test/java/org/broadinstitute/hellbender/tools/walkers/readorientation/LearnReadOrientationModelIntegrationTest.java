package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

public class LearnReadOrientationModelIntegrationTest extends CommandLineProgramTest {
    public static final String testDir = toolsTestDir + "read_orientation_filter/";
    public static final String hapmapBamSnippet = testDir + "hapmap-20-plex-chr-20-21-read-orientation.bam";
    public static final String intervalList = testDir + "hapmap-20-plex-chr-20-21-read-orientation.interval_list";

    @DataProvider(name = "scatterCounts")
    public Object[][] getScatterCounts(){
        return new Object[][]{{1}, {5}};
    }
    /**
     * Test that the tool sites of orientation bias that are manually picked out.
     * Also tests scattering CollectF1R2Counts
     */
    @Test(dataProvider = "scatterCounts")
    public void testOnRealBam(final int scatterCount) throws IOException {
        final File scatteredDir = createTempDir("scattered");

        // Step 1: SplitIntervals
        final File intervalDir = createTempDir("intervals");
        runCommandLine(
            Arrays.asList(
                    "-R", b37_reference_20_21,
                    "-L", intervalList,
                    "-O", intervalDir.getAbsolutePath(),
                    "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount)),
            SplitIntervals.class.getSimpleName());

        // Step 2: CollectF1R2Counts
        final File[] intervals = intervalDir.listFiles();
        final List<File> extractedDirs = IntStream.range(0, intervals.length).mapToObj(i -> createTempDir("extracted_" + i)).collect(Collectors.toList());
        final List<File> scatteredTarGzs = IntStream.range(0, intervals.length).mapToObj(i -> new File(scatteredDir, "scatter_" + i + ".tar.gz")).collect(Collectors.toList());
        for (int i = 0; i < intervals.length; i++){

            runCommandLine(Arrays.asList(
                    "-R", b37_reference_20_21,
                    "-I", hapmapBamSnippet,
                    "-L", intervals[i].getAbsolutePath(),
                    "-O", scatteredTarGzs.get(i).getAbsolutePath()),
                    CollectF1R2Counts.class.getSimpleName());

            IOUtils.extractTarGz(scatteredTarGzs.get(i).toPath(), extractedDirs.get(i).toPath());

            final int iFinal = i;
            final File refHist = F1R2CountsCollector.getRefHistogramsFromExtractedTar(extractedDirs.get(i)).get(0);

            // Ensure that we print every bin, even when the count is 0
            final int lineCount = (int) Files.lines(Paths.get(refHist.getAbsolutePath())).filter(l -> l.matches("^[0-9].+")).count();
            Assert.assertEquals(lineCount, F1R2FilterConstants.DEFAULT_MAX_DEPTH);
        }

        // Step 3: LearnReadOrientationModel
        final File priorTarGz = createTempFile("prior", ".tar.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, priorTarGz.getAbsolutePath());
        IntStream.range(0, intervals.length).forEach(n -> args.addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, scatteredTarGzs.get(n).getAbsolutePath()));

        runCommandLine(args.getArgsList(), LearnReadOrientationModel.class.getSimpleName());

        final File extractedPriorDir = createTempDir("extracted_priors");
        IOUtils.extractTarGz(priorTarGz.toPath(), extractedPriorDir.toPath());


        final ArtifactPriorCollection artifactPriorCollection = ArtifactPriorCollection.readArtifactPriors(extractedPriorDir.listFiles()[0]);

        // Step 4: Mutect 2
        final File unfilteredVcf = GATKBaseTest.createTempFile("unfiltered", ".vcf");
        final File filteredVcf = GATKBaseTest.createTempFile("filtered", ".vcf");
        final File bamout = GATKBaseTest.createTempFile("SM-CEMAH", ".bam");

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                    "-I", hapmapBamSnippet,
                    "-R", b37_reference_20_21,
                    "-O", unfilteredVcf.getAbsolutePath(),
                    "-bamout", bamout.getAbsolutePath()),
                Mutect2.class.getSimpleName()));

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                        "-V", unfilteredVcf.getAbsolutePath(),
                        "-R", b37_reference_20_21,
                        "--" + M2FiltersArgumentCollection.ARTIFACT_PRIOR_TABLE_NAME, priorTarGz.getAbsolutePath(),
                        "-O", filteredVcf.getAbsolutePath()),
                FilterMutectCalls.class.getSimpleName()));


        // These artifacts have been verified manually
        // The pair is of type (Position, Expected Source of Prior Probability)
        // Prior for e.g. TGA->A F2R1 should come from TCA->T F1R2
        final List<Triple<Integer, ReadOrientation, ArtifactState>> knownArtifacts = Arrays.asList(
                new ImmutableTriple<>(23421079, ReadOrientation.F2R1, ArtifactState.F2R1_G), // CAC->G F2R1
                new ImmutableTriple<>(34144749, ReadOrientation.F2R1, ArtifactState.F2R1_A), // TGA->A F2R1, ref context is not canonical
                new ImmutableTriple<>(62165528, ReadOrientation.F1R2, ArtifactState.F1R2_A)); // CGC->A F1R2

        for (final Triple<Integer, ReadOrientation, ArtifactState> artifact : knownArtifacts) {
            final int position = artifact.getLeft();

            Optional<VariantContext> variant = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                    .filter(vc -> vc.getStart() == position).findFirst();
            Assert.assertTrue(variant.isPresent());

            // Check that the expected filters were applied
            Assert.assertTrue(variant.get().getFilters().contains(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME));
        }
    }

    @Test
    public void testTwoSamples() throws Exception {
        final File countsTarGz = createTempFile("counts", ".tar.gz");
        final File priorsTarGz = createTempFile("priors", ".tar.gz");
        final String sample1 = "SAMPLE1";
        final String sample2 = "SAMPLE2";
        final File sam1 = CollectF1R2CountsIntegrationTest.createSyntheticSam(10, 1, sample1);
        final File sam2 = CollectF1R2CountsIntegrationTest.createSyntheticSam(20, 2, sample2);


        new Main().instanceMain(makeCommandLineArgs(Arrays.asList(
                "-R", hg19_chr1_1M_Reference, "-I", sam1.getAbsolutePath(), "-I", sam2.getAbsolutePath(), "-O", countsTarGz.getAbsolutePath()),
                CollectF1R2Counts.class.getSimpleName()));

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, countsTarGz.getAbsolutePath())
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, priorsTarGz.getAbsolutePath());

        runCommandLine(args.getArgsList(), LearnReadOrientationModel.class.getSimpleName());

        final File extractedPriorsDir = createTempDir("extracted");
        IOUtils.extractTarGz(priorsTarGz.toPath(), extractedPriorsDir.toPath());

        Assert.assertTrue(new File(extractedPriorsDir, sample1 + LearnReadOrientationModel.ARTIFACT_PRIOR_EXTENSION).exists());
        Assert.assertTrue(new File(extractedPriorsDir, sample2 + LearnReadOrientationModel.ARTIFACT_PRIOR_EXTENSION).exists());
    }
}