package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;

public class GenotypeGVCFsIntegrationTest extends CommandLineProgramTest {

    //This value was originally 30 but was changed to 10.  We're setting most tests to use 30 to temporarily avoid having
    //to update the expected outputs
    private static final List<String> STAND_CALL_CONF_30 = Arrays.asList("-stand_call_conf", "30");

    private static <T> void assertForEachElementInLists(final List<T> actual, final List<T> expected, final BiConsumer<T, T> assertion) {
        Assert.assertEquals(actual.size(), expected.size(), "different number of elements in lists");
        for (int i = 0; i < actual.size(); i++) {
            assertion.accept(actual.get(i), expected.get(i));
        }
    }

    @DataProvider(name = "gvcfsToGenotype")
    public Object[][] gvcfsToGenotype() {
        String basePairGVCF = "gvcf.basepairResolution.gvcf";
        return new Object[][]{
                {basePairGVCF, "gvcf.basepairResolution.output.vcf", STAND_CALL_CONF_30}, //base pair level gvcf
                {"testUpdatePGT.gvcf", "testUpdatePGT.output.vcf", STAND_CALL_CONF_30},   //testUpdatePGTStrandAlleleCountsBySample
                {"gvcfExample1.vcf", "gvcfExample1.vcf.expected.vcf", STAND_CALL_CONF_30}, //single sample vcf
                {"combined_genotype_gvcf_exception.vcf", "combined_genotype_gvcf_exception.output.vcf", STAND_CALL_CONF_30}, //test that an input vcf with 0/0 already in GT field is overwritten
                {"combined_genotype_gvcf_exception.nocall.vcf", "combined_genotype_gvcf_exception.output.vcf", STAND_CALL_CONF_30},  //same test as above but with ./.
                {basePairGVCF, "ndaTest.expected.vcf", Arrays.asList("-nda", "-stand_call_conf", "30")},  //annotating with the number of alleles discovered option
                {basePairGVCF, "maxAltAllelesTest.expected.vcf", Arrays.asList("--maxAltAlleles", "1", "-stand_call_conf", "30")}, //restricting the max number of alt alleles
                {basePairGVCF, "standardConfTest.expected.vcf", Arrays.asList("-stand_call_conf", "300")}, //changing call confidence
                {"spanningDel.combined.g.vcf", "spanningDel.combined.g.vcf.expected.vcf", STAND_CALL_CONF_30},
                {"spanningDel.delOnly.g.vcf", "spanningDel.delOnly.g.vcf.expected.vcf", STAND_CALL_CONF_30},
                {"CEUTrio.20.21.gatk3.4.g.vcf", "CEUTrio.20.21.expected.vcf", Arrays.asList("--dbsnp", "src/test/resources/large/dbsnp_138.b37.20.21.vcf", "-stand_call_conf", "30")},
                //{basePairGVCF, "gvcf.basepairResolution.includeNonVariantSites.expected.vcf", Collections.singletonList("--includeNonVariantSites")
        };
    }

    @Test(dataProvider = "gvcfsToGenotype")
    public void testGenotypesOnly(String input, String expected, List<String> extraArgs) throws IOException {
        assertGenotypesMatch(getTestFile(input), getTestFile(expected), extraArgs);
    }

    @Test(dataProvider = "gvcfsToGenotype")
    public void testEntireVariantContext(String input, String expected, List<String> extraArgs) throws IOException {
        assertVariantContextsMatch(getTestFile(input), getTestFile(expected), extraArgs);
    }

    private void assertVariantContextsMatch(File input, File expected, List<String> extraArgs) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input, expected, extraArgs, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                Arrays.asList("FS", //TODO There's a bug in GATK 3 computing FS
                        "QD", //TODO QD has a cap value and anything that reaches that is randomized.  It's difficult to reproduce the same random numbers accross gatk3 -> 4
                        "InbreedingCoeff")));
    }

    private void assertGenotypesMatch(File input, File expected, List<String> additionalArguments) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input, expected, additionalArguments, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes);
    }

    private void runGenotypeGVCFSAndAssertSomething(File input, File expected, List<String> additionalArguments, BiConsumer<VariantContext, VariantContext> assertion) throws IOException {
        final File output = createTempFile("genotypegvcf", "vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addVCF(input)
                .addOutput(output);

        additionalArguments.forEach(args::add);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, assertion); //TODO Inbreeding calculation changed between 3.4 and now
    }

    /**
     * Returns a list of VariantContext records from a VCF file
     *
     * @param vcfFile VCF file
     * @return list of VariantContext records
     * @throws IOException if the file does not exist or can not be opened
     */
    private static List<VariantContext> getVariantContexts(final File vcfFile) throws IOException {
        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(vcfFile);
        final LineIterator lineIteratorVCF = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIteratorVCF);

        final List<VariantContext> VCs = new ArrayList<>();
        while (lineIteratorVCF.hasNext()) {
            final String line = lineIteratorVCF.next();
            Assert.assertFalse(line == null);
            VCs.add(codec.decode(line));
        }

        return VCs;
    }
}
