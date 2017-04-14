package org.broadinstitute.hellbender.engine;

import com.intel.genomicsdb.GenomicsDBUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;


public class GenomicsDBIntegrationTest extends CommandLineProgramTest {
    public static final String TEST_DATA_PATH =  publicTestDir + "/org/broadinstitute/hellbender/engine/GenomicsDBIntegration/";
    public static final File LOADER_JSON = new File(TEST_DATA_PATH, "tiny_jsons/loader.json");

    private static final String TINY_GVCF = TEST_DATA_PATH + "tiny_jsons/tiny.gvcf";

    private static final String WORKSPACE = TEST_DATA_PATH + "workspace";
    private static final String GENOMICS_DB_JSONS = "gendb://" + TEST_DATA_PATH + "tiny_jsons";
    
    @Override
    public String getTestedClassName() {
        return SelectVariants.class.getSimpleName();
    }

    @Test
    public void testGenomicsDBInClassPath(){
        final String path = "/genomicsdb/"+System.mapLibraryName("tiledbgenomicsdb");
        Assert.assertNotNull(GenomicsDBUtils.class.getResource(path), "Could not find the genomicsdb binary at " + path);
    }

    @Test
    public void testAsDrivingVariants() throws IOException {
        final File expected = new File(TINY_GVCF);
        testExpectedVariantsFromGenomicsDB(expected, new ArgumentsBuilder()
                .addArgument("V", GENOMICS_DB_JSONS));
    }

    @Test
    public void testWithIntervalsAsAuxiliary() throws IOException {
        final File expected = new File(TINY_GVCF);
        testExpectedVariantsFromGenomicsDB(expected, new ArgumentsBuilder()
                .addArgument("V", TINY_GVCF)
            .addArgument("concordance", GENOMICS_DB_JSONS)
            .addArgument("L", "20")
            .addReference(new File(largeFileTestDir,"human_g1k_v37.20.21.fasta")));
    }

    @Test
    public void testWithIntervalsDriving() throws IOException {
        final File expected = new File(TINY_GVCF);
        testExpectedVariantsFromGenomicsDB(expected, new ArgumentsBuilder()
                .addArgument("V", GENOMICS_DB_JSONS)
                .addArgument("L", "20")
                .addReference(new File(largeFileTestDir, "human_g1k_v37.20.21.fasta")));
    }

    @Test
    public void testRestrictingIntervals() throws IOException {
        testExpectedVariantsFromGenomicsDB(new File(TEST_DATA_PATH, "intervalsRestrictedExpected.gvcf"), new ArgumentsBuilder()
                .addArgument("V", GENOMICS_DB_JSONS)
                .addReference(new File(largeFileTestDir, "human_g1k_v37.20.21.fasta"))
                .addArgument("L", "20:69491-69521"));
    }

    private void testExpectedVariantsFromGenomicsDB(File expected, ArgumentsBuilder baseArgs) throws IOException {
        final File output = createTempFile("variants", ".vcf");
        GenomicsDBTestUtils.runOnGenomicsDBArray( () -> {
            final ArgumentsBuilder args = baseArgs
                    .addOutput(output);
            runCommandLine(args);
        }, new File(WORKSPACE), LOADER_JSON);

        final Iterable<VariantContext> actualVcs = new FeatureDataSource<>(output);
        final Iterable<VariantContext> expectedVcs = new FeatureDataSource<>(expected);

        GenomicsDBTestUtils.assertCondition(actualVcs, expectedVcs, (a,e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a,e, Collections.emptyList()));
    }

    @Test
    public void testAsAuxiliaryVariants() throws IOException {
        testExpectedVariantsFromGenomicsDB(new File(TINY_GVCF), new ArgumentsBuilder()
                    .addArgument("V", TINY_GVCF)
                    .addArgument("concordance", GENOMICS_DB_JSONS));
    }
}
