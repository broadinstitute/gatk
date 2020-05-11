package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3BaseData;
import org.broadinstitute.hellbender.CommandLineProgramTest;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;


public class GeneExpressionEvaluationIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert a match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=GeneExpressionEvaluationIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_MATCH_EXPECTED_OUTPUTS = false;

    public static final String TEST_FILES_DIR = toolsTestDir + "walkers/rnaseq/GeneExpressionEvaluation/";

    //allow results to differ by less than 0.1% from past and still count as passing
    private static final double DEFAULT_LENIENCE = 0.001;

    /*
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @DataProvider(name = "consistentWithPastDataProvider")
    public Object[][] consistentWithPastDataProvider() {
        return new Object[][] {
                {GeneExpressionEvaluation.MultiOverlapMethod.PROPORTIONAL, GeneExpressionEvaluation.MultiMapMethod.IGNORE, new File(TEST_FILES_DIR, "expected.overlap.proportional.map.ignore.tsv")},
                {GeneExpressionEvaluation.MultiOverlapMethod.EQUAL, GeneExpressionEvaluation.MultiMapMethod.IGNORE, new File(TEST_FILES_DIR, "expected.overlap.equal.map.ignore.tsv")},
                {GeneExpressionEvaluation.MultiOverlapMethod.PROPORTIONAL, GeneExpressionEvaluation.MultiMapMethod.EQUAL, new File(TEST_FILES_DIR, "expected.overlap.proportional.map.equal.tsv")},
                {GeneExpressionEvaluation.MultiOverlapMethod.EQUAL, GeneExpressionEvaluation.MultiMapMethod.EQUAL, new File(TEST_FILES_DIR, "expected.overlap.equal.map.equal.tsv")}
        };
    }

    @Test(dataProvider = "consistentWithPastDataProvider")
    public void testConsistentWithPast(final GeneExpressionEvaluation.MultiOverlapMethod multiOverlapMethod, final GeneExpressionEvaluation.MultiMapMethod multiMapMethod, final File expected) throws IOException {
        final File output = createTempFile("testConsistentWithPast", ".tsv");

        final String outputPath = UPDATE_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", NA12878_20_RNAseq_bam,
                "-N", "NA12878",
                "-G", b37_20_gff3,
                "--gene_id_key", "Name",
                "--grouping_type", "gene",
                "--grouping_type", "pseudogene",
                "--overlap_type", "exon",
                "--multiOverlapMethod", multiOverlapMethod.toString(),
                "--multiMapMethod", multiMapMethod.toString(),
                "-O", outputPath
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_MATCH_EXPECTED_OUTPUTS ) {
            assertResultsEquivalent(output.toPath(), expected.toPath(), DEFAULT_LENIENCE);
        }
    }

    @Test
    public void testSwapTranscriptionRead() throws IOException {
        final File output1 = createTempFile("testSwapTranscriptionRead1", ".tsv");
        final File output2 = createTempFile("testSwapTranscriptionRead2", ".tsv");

        final String[] args1 = {
                "-I", NA12878_20_RNAseq_bam,
                "-N", "NA12878",
                "-G", b37_20_gff3,
                "--gene_id_key", "Name",
                "--grouping_type", "gene",
                "--grouping_type", "pseudogene",
                "--overlap_type", "exon",
                "--trancriptionRead", "R1",
                "-O", output1.getAbsolutePath()
        };

        final String[] args2 = {
                "-I", NA12878_20_RNAseq_bam,
                "-N", "NA12878",
                "-G", b37_20_gff3,
                "--gene_id_key", "Name",
                "--grouping_type", "gene",
                "--grouping_type", "pseudogene",
                "--overlap_type", "exon",
                "--trancriptionRead", "R2",
                "-O", output2.getAbsolutePath()
        };

        runCommandLine(args1);
        runCommandLine(args2);
        assertResultsEquivalent(output1.toPath(), output2.toPath(), 0, true);
    }

    private void assertResultsEquivalent(final Path file1, final Path file2, final double lenience) throws IOException {
        assertResultsEquivalent(file1, file2, lenience, false);
    }

    private void assertResultsEquivalent(final Path file1, final Path file2, final double lenience, final boolean swapTranscriptionStrand) throws IOException {
        final Map<Gff3BaseData, GeneExpressionEvaluation.Coverage> coverageMap1 = readCoverage(file1);
        final Map<Gff3BaseData, GeneExpressionEvaluation.Coverage> coverageMap2 = readCoverage(file2);

        Assert.assertEquals(coverageMap1.size(), coverageMap2.size());

        for (final Map.Entry<Gff3BaseData, GeneExpressionEvaluation.Coverage> entry : coverageMap1.entrySet()) {
            final Gff3BaseData feature = entry.getKey();
            final GeneExpressionEvaluation.Coverage coverage1 = entry.getValue();

            Assert.assertTrue(coverageMap2.containsKey(feature));
            final GeneExpressionEvaluation.Coverage coverage2 = coverageMap2.get(feature);

            if (feature.getStrand() != Strand.NONE) {
                assertEquivalentWithLenience(coverage1.getSenseCount(), swapTranscriptionStrand? coverage2.getAntisenseCount() : coverage2.getSenseCount(), lenience);
                assertEquivalentWithLenience(coverage2.getSenseCount(), swapTranscriptionStrand? coverage1.getAntisenseCount() : coverage1.getSenseCount(), lenience);
            } else {
                assertEquivalentWithLenience(coverage1.getSenseCount(), coverage2.getSenseCount(), lenience);
                Assert.assertEquals(coverage1.getAntisenseCount(), 0);
                Assert.assertEquals(coverage2.getAntisenseCount(), 0);
            }
        }
    }

    private void assertEquivalentWithLenience(final float res1, final float res2, final double lenience) {
        Assert.assertTrue(res1 >= 0);
        Assert.assertTrue(res2 >= 0);
        if (res1 + res2 > 0) {
            Assert.assertTrue(2*Math.abs(res1 - res2)/(res1 + res2) <= lenience, "res1 = " + res1 + " res2 = " + res2);
        }
    }

    private Map<Gff3BaseData, GeneExpressionEvaluation.Coverage> readCoverage(final Path file) throws IOException {
        final GeneExpressionEvaluation.FragmentCountReader reader = new GeneExpressionEvaluation.FragmentCountReader(file);
        final Map<Gff3BaseData, GeneExpressionEvaluation.Coverage> coverageMap = new LinkedHashMap<>();

        for (final GeneExpressionEvaluation.SingleStrandFeatureCoverage singleStrandFeatureCoverage : reader) {
            coverageMap.compute(
                    singleStrandFeatureCoverage.baseData, (k, v) -> v == null ?
                            (singleStrandFeatureCoverage.sense? new GeneExpressionEvaluation.Coverage(singleStrandFeatureCoverage.count, 0) :
                                    new GeneExpressionEvaluation.Coverage(0, singleStrandFeatureCoverage.count)
                            ) :
                            (singleStrandFeatureCoverage.sense? new GeneExpressionEvaluation.Coverage(singleStrandFeatureCoverage.count, v.getAntisenseCount()) :
                                    new GeneExpressionEvaluation.Coverage(v.getSenseCount(), singleStrandFeatureCoverage.count)
                            )
            );
        }
        return coverageMap;
    }
}

