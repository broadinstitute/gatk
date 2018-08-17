package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.ClinVarFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.LmmFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.LofFilter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class FilterFuncotationsIntegrationTest extends CommandLineProgramTest {

    private static final Path TEST_DATA_DIR = getTestDataDir().toPath().resolve("FilterFuncotations");

    private static final Set<String> ALL_FILTERS = new HashSet<>(Arrays.asList(
            ClinVarFilter.CLINSIG_INFO_VALUE, LofFilter.CLINSIG_INFO_VALUE, LmmFilter.CLINSIG_INFO_VALUE));

    @DataProvider(name = "uniformVcfProvider")
    public Object[][] uniformVcfProvider() {
        return new Object[][]{
                {"clinvar.vcf", FilterFuncotations.Reference.hg19, Collections.emptySet(), Collections.singleton(ClinVarFilter.CLINSIG_INFO_VALUE)},
                {"lmm.vcf", FilterFuncotations.Reference.hg38, Collections.emptySet(), Collections.singleton(LmmFilter.CLINSIG_INFO_VALUE)},
                {"lof.vcf", FilterFuncotations.Reference.b37, Collections.emptySet(), Collections.singleton(LofFilter.CLINSIG_INFO_VALUE)},
                {"all.vcf", FilterFuncotations.Reference.hg38, Collections.emptySet(), ALL_FILTERS},
                {"multi-transcript.vcf", FilterFuncotations.Reference.hg38, Collections.emptySet(), ALL_FILTERS},
                {"multi-allelic.vcf", FilterFuncotations.Reference.hg38, Collections.emptySet(), ALL_FILTERS},
                {"none.vcf", FilterFuncotations.Reference.hg38, Collections.singleton(FilterFuncotationsConstants.NOT_CLINSIG_FILTER),
                        Collections.singleton(FilterFuncotationsConstants.CLINSIG_INFO_NOT_SIGNIFICANT)}
        };
    }

    @Test(dataProvider = "uniformVcfProvider")
    public void testFilterUniform(final String vcfName,
                                  final FilterFuncotations.Reference ref,
                                  final Set<String> expectedFilters,
                                  final Set<String> expectedAnnotations) {

        final File tmpOut = createTempFile(vcfName + ".filtered", ".vcf");

        final List<String> args = Arrays.asList(
                "-V", TEST_DATA_DIR.resolve(vcfName).toString(),
                "-O", tmpOut.toString(),
                "--ref-version", ref.name()
        );
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(tmpOut.toString());
        vcf.getRight().forEach(variant -> {
            Assert.assertEquals(variant.getFilters(), expectedFilters);

            final List<String> clinsigAnnotations = variant.getCommonInfo()
                    .getAttributeAsStringList(FilterFuncotationsConstants.CLINSIG_INFO_KEY, "");
            Assert.assertEquals(new HashSet<>(clinsigAnnotations), expectedAnnotations);
        });
    }
}
