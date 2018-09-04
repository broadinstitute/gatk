package org.broadinstitute.hellbender.engine;

import com.intel.genomicsdb.GenomicsDBLibLoader;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;


public class GenomicsDBIntegrationTest extends CommandLineProgramTest {
    public static final String TEST_DATA_PATH =  publicTestDir + "/org/broadinstitute/hellbender/engine/GenomicsDBIntegration/";
    private static final File TINY_GVCF = new File(TEST_DATA_PATH, "tiny.g.vcf");
    private static final SimpleInterval INTERVAL = new SimpleInterval("20", 1, 63025520);

    @Override
    public String getTestedClassName() {
        return SelectVariants.class.getSimpleName();
    }

    @Test
    public void testGenomicsDBInClassPath(){
        final String path = "/"+System.mapLibraryName("tiledbgenomicsdb");
        Assert.assertNotNull(GenomicsDBLibLoader.class.getResource(path), "Could not find the genomicsdb binary at " + path);
    }

    @Test
    public void testAsDrivingVariants() throws IOException {
        final File workspace = GenomicsDBTestUtils.createTempGenomicsDB(TINY_GVCF, INTERVAL);
        testExpectedVariantsFromGenomicsDB(TINY_GVCF, new ArgumentsBuilder()
                .addArgument("V", GenomicsDBTestUtils.makeGenomicsDBUri(workspace))
                .addReference(new File(GATKBaseTest.b37_reference_20_21)));
    }

    @Test
    public void testWithIntervalsAsAuxiliary() throws IOException {
        final File workspace = GenomicsDBTestUtils.createTempGenomicsDB(TINY_GVCF, INTERVAL);

        testExpectedVariantsFromGenomicsDB(TINY_GVCF, new ArgumentsBuilder()
                .addArgument("V", TINY_GVCF.getAbsolutePath())
            .addArgument("concordance", GenomicsDBTestUtils.makeGenomicsDBUri(workspace))
            .addArgument("L", "20")
            .addReference(new File(GATKBaseTest.b37_reference_20_21)));
    }

    @Test
    public void testWithIntervalsDriving() throws IOException {
        final File workspace = GenomicsDBTestUtils.createTempGenomicsDB(TINY_GVCF, INTERVAL);
        testExpectedVariantsFromGenomicsDB(TINY_GVCF, new ArgumentsBuilder()
                .addArgument("L", "20")
                .addReference(new File(GATKBaseTest.b37_reference_20_21))
                .addArgument("V", GenomicsDBTestUtils.makeGenomicsDBUri(workspace))
        );
    }


    @Test
    public void testRestrictingIntervals() throws IOException {
        final File workspace = GenomicsDBTestUtils.createTempGenomicsDB(TINY_GVCF, INTERVAL);
        testExpectedVariantsFromGenomicsDB(new File(TEST_DATA_PATH, "intervalsRestrictedExpected.g.vcf"), new ArgumentsBuilder()
                .addReference(new File(GATKBaseTest.b37_reference_20_21))
                .addArgument("L", "20:69491-69521")
                .addArgument("V", GenomicsDBTestUtils.makeGenomicsDBUri(workspace)));
    }

    private void testExpectedVariantsFromGenomicsDB(File expected, ArgumentsBuilder baseArgs) throws IOException {
        final File output = createTempFile("variants", ".vcf");
        final ArgumentsBuilder args = baseArgs
                .addOutput(output);
        runCommandLine(args);

        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(expected)) {
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                                     (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                                                                                                     Collections.emptyList()));
        }
    }

    @Test
    public void testAsAuxiliaryVariants() throws IOException {
        final File workspace = GenomicsDBTestUtils.createTempGenomicsDB(TINY_GVCF, INTERVAL);
        testExpectedVariantsFromGenomicsDB(TINY_GVCF, new ArgumentsBuilder()
                    .addArgument("V", TINY_GVCF.getAbsolutePath())
                    .addArgument("concordance", GenomicsDBTestUtils.makeGenomicsDBUri(workspace))
                    .addReference(new File(GATKBaseTest.b37_reference_20_21)));
    }
}
