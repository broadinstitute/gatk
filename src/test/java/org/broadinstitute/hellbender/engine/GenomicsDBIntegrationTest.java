package org.broadinstitute.hellbender.engine;

import com.intel.genomicsdb.GenomicsDBUtils;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;


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
        Assert.assertNotNull(GenomicsDBUtils.class.getResource(path), "Could not find the genomicsdb binary at " + path);
    }

    @Test
    public void testAsDrivingVariants() throws IOException {
        final String workspace = createTempGenomicsDB(TINY_GVCF, INTERVAL);
        testExpectedVariantsFromGenomicsDB(TINY_GVCF, new ArgumentsBuilder()
                .addArgument("V", makeGenDBUri(workspace))
                .addReference(new File(BaseTest.b37_reference_20_21)));
    }

    @Test
    public void testWithIntervalsAsAuxiliary() throws IOException {
        final String workspace = createTempGenomicsDB(TINY_GVCF,INTERVAL);

        testExpectedVariantsFromGenomicsDB(TINY_GVCF, new ArgumentsBuilder()
                .addArgument("V", TINY_GVCF.getAbsolutePath())
            .addArgument("concordance", makeGenDBUri(workspace))
            .addArgument("L", "20")
            .addReference(new File(BaseTest.b37_reference_20_21)));
    }

    @Test
    public void testWithIntervalsDriving() throws IOException {
        final String workspace = createTempGenomicsDB(TINY_GVCF, INTERVAL);
        testExpectedVariantsFromGenomicsDB(TINY_GVCF, new ArgumentsBuilder()
                .addArgument("L", "20")
                .addReference(new File(BaseTest.b37_reference_20_21))
                .addArgument("V", makeGenDBUri(workspace))
        );
    }


    private static String makeGenDBUri(String workspace){
        return FeatureDataSource.GENOMIC_DB_URI_SCHEME + workspace;
    }

    @Test
    public void testRestrictingIntervals() throws IOException {
        final String workspace = createTempGenomicsDB(TINY_GVCF, INTERVAL);
        testExpectedVariantsFromGenomicsDB(new File(TEST_DATA_PATH, "intervalsRestrictedExpected.g.vcf"), new ArgumentsBuilder()
                .addReference(new File(BaseTest.b37_reference_20_21))
                .addArgument("L", "20:69491-69521")
                .addArgument("V", makeGenDBUri(workspace)));
    }

    private void testExpectedVariantsFromGenomicsDB(File expected, ArgumentsBuilder baseArgs) throws IOException {
        final File output = createTempFile("variants", ".vcf");
        final ArgumentsBuilder args = baseArgs
                .addOutput(output);
        runCommandLine(args);

        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(expected)) {
            BaseTest.assertCondition(actualVcs, expectedVcs,
                                     (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                                                                                                     Collections.emptyList()));
        }
    }

    private static String createTempGenomicsDB(final File gvcfs, final Locatable interval) {
        final File workspaceDir = BaseTest.createTempDir("genomicsDBWorkspace");

        final CommandLineProgramTest importer = new CommandLineProgramTest() {
            @Override
            public String getTestedToolName() {
                return GenomicsDBImport.class.getSimpleName();
            }
        };

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(gvcfs);

        final String workspace = new File(workspaceDir, "workspace").getAbsolutePath();
        args.addArgument(GenomicsDBImport.WORKSPACE_ARG_NAME, workspace);
        args.addArgument("L", IntervalUtils.locatableToString(INTERVAL));
        importer.runCommandLine(args);
        return workspace;
    }

    @Test
    public void testAsAuxiliaryVariants() throws IOException {
        final String workspace = createTempGenomicsDB(TINY_GVCF, INTERVAL);
        testExpectedVariantsFromGenomicsDB(TINY_GVCF, new ArgumentsBuilder()
                    .addArgument("V", TINY_GVCF.getAbsolutePath())
                    .addArgument("concordance", makeGenDBUri(workspace))
                    .addReference(new File(BaseTest.b37_reference_20_21)));
    }
}
