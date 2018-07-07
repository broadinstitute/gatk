package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class FilterFuncotationsIntegrationTest extends CommandLineProgramTest {

    private static final Path TEST_DATA_DIR = getTestDataDir().toPath().resolve("FilterFuncotations");
    
    @DataProvider(name = "uniformVcfProvider")
    public Object[][] uniformVcfProvider() {
        return new Object[][]{
                {"clinvar.vcf", 19, Collections.emptySet(), Collections.singleton("CLINVAR")},
                {"lmm.vcf", 38, Collections.emptySet(), Collections.singleton("LMM")},
                {"lof.vcf", 19, Collections.emptySet(), Collections.singleton("LOF")},
                {"all.vcf", 38, Collections.emptySet(), new HashSet<>(Arrays.asList("CLINVAR", "LMM", "LOF"))},
                {"none.vcf", 38, Collections.singleton(FilterFuncotations.NOT_CLINSIG_FILTER), Collections.singleton("NONE")}
        };
    }

    @Test(dataProvider = "uniformVcfProvider")
    public void testFilterUniform(final String vcfName,
                                  final int build,
                                  final Set<String> expectedFilters,
                                  final Set<String> expectedAnnotations) throws IOException {

        final Path tmpOut = Files.createTempFile(vcfName + ".filtered", ".vcf");
        IOUtil.deleteOnExit(tmpOut);

        final List<String> args = Arrays.asList(
                "-V", TEST_DATA_DIR.resolve(vcfName).toString(),
                "-O", tmpOut.toString(),
                "--ref-version", "hg" + build
        );
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(tmpOut.toString());
        vcf.getRight().forEach(variant -> {
            Assert.assertEquals(variant.getFilters(), expectedFilters);

            final List<String> clinsigAnnotations = variant.getCommonInfo()
                    .getAttributeAsStringList(FilterFuncotations.CLINSIG_RULE_KEY, "");
            Assert.assertEquals(new HashSet<>(clinsigAnnotations), expectedAnnotations);
        });
    }
}
