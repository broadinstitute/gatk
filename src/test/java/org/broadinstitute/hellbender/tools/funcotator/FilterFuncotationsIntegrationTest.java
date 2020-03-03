package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.AutosomalRecessiveConstants;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.ClinVarFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.LmmFilter;
import org.broadinstitute.hellbender.tools.funcotator.filtrationRules.LofFilter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.util.*;

public class FilterFuncotationsIntegrationTest extends CommandLineProgramTest {

    private static final Path TEST_DATA_DIR = getTestDataDir().toPath().resolve("FilterFuncotations");

    private static final Set<String> ALL_FILTERS = new HashSet<>(Arrays.asList(
            ClinVarFilter.CLINSIG_INFO_VALUE, LofFilter.CLINSIG_INFO_VALUE, LmmFilter.CLINSIG_INFO_VALUE, AutosomalRecessiveConstants.AR_INFO_VALUE));

    @DataProvider(name = "uniformVcfProvider")
    public Object[][] uniformVcfProvider() {
        return new Object[][]{
                {"all.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.emptySet(), ALL_FILTERS},
                {"all_gnomad.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.gnomad, Collections.emptySet(), ALL_FILTERS},
                {"ar_homvar.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.emptySet(), Collections.singleton(AutosomalRecessiveConstants.AR_INFO_VALUE)},
                {"ar_hetvar.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.singleton(FilterFuncotationsConstants.NOT_CLINSIG_FILTER), Collections.singleton(FilterFuncotationsConstants.CLINSIG_INFO_NOT_SIGNIFICANT)},
                {"ar_compound_het.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.emptySet(), Collections.singleton(AutosomalRecessiveConstants.AR_INFO_VALUE)},
                {"clinvar.vcf", FilterFuncotations.Reference.hg19, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.emptySet(), Collections.singleton(ClinVarFilter.CLINSIG_INFO_VALUE)},
                {"clinvar_gnomad.vcf", FilterFuncotations.Reference.hg19, FilterFuncotations.AlleleFrequencyDataSource.gnomad, Collections.emptySet(), Collections.singleton(ClinVarFilter.CLINSIG_INFO_VALUE)},
                {"gnomad_af_failing_cases.vcf", FilterFuncotations.Reference.hg19, FilterFuncotations.AlleleFrequencyDataSource.gnomad, Collections.singleton(FilterFuncotationsConstants.NOT_CLINSIG_FILTER),
                        Collections.singleton(FilterFuncotationsConstants.CLINSIG_INFO_NOT_SIGNIFICANT)},
                {"gnomad_af_passing_cases.vcf", FilterFuncotations.Reference.hg19, FilterFuncotations.AlleleFrequencyDataSource.gnomad, Collections.emptySet(), Collections.singleton(LofFilter.CLINSIG_INFO_VALUE)},
                {"lmm.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.emptySet(), Collections.singleton(LmmFilter.CLINSIG_INFO_VALUE)},
                {"lof.vcf", FilterFuncotations.Reference.b37, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.emptySet(), Collections.singleton(LofFilter.CLINSIG_INFO_VALUE)},
                {"lof_gnomad.vcf", FilterFuncotations.Reference.b37, FilterFuncotations.AlleleFrequencyDataSource.gnomad, Collections.emptySet(), Collections.singleton(LofFilter.CLINSIG_INFO_VALUE)},
                {"multi-allelic.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.emptySet(), ALL_FILTERS},
                {"multi-allelic_gnomad.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.gnomad, Collections.emptySet(), ALL_FILTERS},
                {"multi-transcript.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.emptySet(), ALL_FILTERS},
                {"multi-transcript_gnomad.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.gnomad, Collections.emptySet(), ALL_FILTERS},
                {"none.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.exac, Collections.singleton(FilterFuncotationsConstants.NOT_CLINSIG_FILTER),
                        Collections.singleton(FilterFuncotationsConstants.CLINSIG_INFO_NOT_SIGNIFICANT)},
                {"none_gnomad.vcf", FilterFuncotations.Reference.hg38, FilterFuncotations.AlleleFrequencyDataSource.gnomad, Collections.singleton(FilterFuncotationsConstants.NOT_CLINSIG_FILTER),
                        Collections.singleton(FilterFuncotationsConstants.CLINSIG_INFO_NOT_SIGNIFICANT)}
        };
    }

    @Test(dataProvider = "uniformVcfProvider")
    public void testFilterUniform(final String vcfName,
                                  final FilterFuncotations.Reference ref,
                                  final FilterFuncotations.AlleleFrequencyDataSource afDataSource,
                                  final Set<String> expectedFilters,
                                  final Set<String> expectedAnnotations) {

        final File tmpOut = createTempFile(vcfName + ".filtered", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add(StandardArgumentDefinitions.VARIANT_SHORT_NAME, TEST_DATA_DIR.resolve(vcfName).toFile())
                .add("ref-version", ref.name())
                .add("allele-frequency-data-source", afDataSource.name())
                .addOutput(tmpOut);

        runCommandLine(args.getArgsArray());

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(tmpOut.toString());
        vcf.getRight().forEach(variant -> {
            Assert.assertEquals(variant.getFilters(), expectedFilters);

            final List<String> clinsigAnnotations = variant.getCommonInfo()
                    .getAttributeAsStringList(FilterFuncotationsConstants.CLINSIG_INFO_KEY, "");
            Assert.assertEquals(new HashSet<>(clinsigAnnotations), expectedAnnotations);
        });
    }
}
