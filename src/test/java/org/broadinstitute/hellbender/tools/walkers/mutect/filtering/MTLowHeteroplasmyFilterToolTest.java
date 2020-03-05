package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class MTLowHeteroplasmyFilterToolTest extends CommandLineProgramTest {
    private static final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");
    private static final File NA12878_MITO_FILTERED_VCF = new File(toolsTestDir, "mutect/mito/filtered.vcf");

    @Test
    public void testLowHetVariantsFiltered() {
        final Set<Integer> low_het_sites = new HashSet<>(Arrays.asList(301, 302));
        final File outputFile = createTempFile("low-het-test", ".vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(MITO_REF.getAbsolutePath())
                .add(StandardArgumentDefinitions.VARIANT_SHORT_NAME, NA12878_MITO_FILTERED_VCF.getAbsolutePath())
                .add(MTLowHeteroplasmyFilterTool.MAX_ALLOWED_LOW_HETS_LONG_NAME, 0)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        Set<VariantContext> variants = VariantContextTestUtils.streamVcf(outputFile)
                .filter(vcf -> vcf.getFilters().contains(GATKVCFConstants.LOW_HET_FILTER_NAME)).collect(Collectors.toSet());
        Assert.assertEquals(variants.stream().map(var -> var.getStart()).collect(Collectors.toList()), low_het_sites, "exprected these sites to have " + GATKVCFConstants.LOW_HET_FILTER_NAME + " filter.");
    }

    @Test
    public void testNoLowHetVariantsFiltered() {
        final File outputFile = createTempFile("no-low-het-test", ".vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addReference(MITO_REF.getAbsolutePath())
                .add(StandardArgumentDefinitions.VARIANT_SHORT_NAME, NA12878_MITO_FILTERED_VCF.getAbsolutePath())
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        Assert.assertTrue(VariantContextTestUtils.streamVcf(outputFile)
                .map(VariantContext::getFilters).noneMatch(filterSet -> filterSet.contains(GATKVCFConstants.LOW_HET_FILTER_NAME)),
                "exprected no variants to have " + GATKVCFConstants.LOW_HET_FILTER_NAME + " filter.");
    }

}
