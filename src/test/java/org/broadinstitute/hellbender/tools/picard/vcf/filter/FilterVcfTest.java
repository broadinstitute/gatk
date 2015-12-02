package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.ListMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Tests for VCF filtration
 */
public final class FilterVcfTest extends CommandLineProgramTest {
    private final File INPUT = new File(getTestDataDir(), "picard/vcf/testFiltering.vcf");

    /** Tests that all records get PASS set as their filter when extreme values are used for filtering. */
    @Test
    public void testNoFiltering() throws Exception {
        final File out = testFiltering(INPUT, 0, 0, 0, Double.MAX_VALUE);
        final VCFFileReader in = new VCFFileReader(out, false);
        for (final VariantContext ctx : in) {
            if (!ctx.filtersWereApplied() || ctx.isFiltered()) {
                Assert.fail("Context should not have been filtered: " + ctx.toString());
            }
        }
    }

    /** Tests that sites with a het allele balance < 0.4 are marked as filtered out. */
    @Test
    public void testAbFiltering() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("tf2", "rs28566954", "rs28548431");
        final File out = testFiltering(INPUT, 0.4, 0, 0, Double.MAX_VALUE);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(filters.keySet(), fails, "Failed sites did not match expected set of failed sites.");
    }

    /** Tests that genotypes with DP < 18 are marked as failed, but not >= 18. */
    @Test
    public void testDpFiltering() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("rs71509448", "rs71628926", "rs13302979", "rs2710876");
        final File out = testFiltering(INPUT, 0, 18, 0, Double.MAX_VALUE);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(filters.keySet(), fails, "Failed sites did not match expected set of failed sites.");
    }

    /** Tests that genotypes with low GQ are filtered appropriately. */
    @Test
    public void testGqFiltering() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("rs71509448"); // SNP with GQ=21; lowest GQ in file

        {
            final File out = testFiltering(INPUT, 0, 0, 20, Double.MAX_VALUE);
            final ListMap<String, String> filters = slurpFilters(out);
            Assert.assertEquals(filters.size(), 0, "Should not have filtered sites: " + filters);
        }
        {
            final File out = testFiltering(INPUT, 0, 0, 21, Double.MAX_VALUE);
            final ListMap<String, String> filters = slurpFilters(out);
            Assert.assertEquals(filters.size(), 0, "Should not have filtered sites: " + filters);
        }
        {
            final File out = testFiltering(INPUT, 0, 0, 22, Double.MAX_VALUE);
            final ListMap<String, String> filters = slurpFilters(out);
            Assert.assertEquals(filters.keySet(), fails, "Failed sites did not match expected set of failed sites.");
        }
    }

    /** Tests that genotypes with DP < 18 are marked as failed, but not >= 18. */
    @Test
    public void testFsFiltering() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("rs13303033", "rs28548431", "rs2799066");
        final File out = testFiltering(INPUT, 0, 0, 0, 5.0d);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(filters.keySet(), fails, "Failed sites did not match expected set of failed sites.");
    }

    @Test
    public void testCombinedFiltering() throws Exception {
        final TreeSet<String> fails = new TreeSet<>(CollectionUtil.makeSet("rs13302979", "rs13303033", "rs2710876", "rs2799066", "rs28548431", "rs28566954", "rs71509448", "rs71628926", "tf2"));
        final File out = testFiltering(INPUT, 0.4, 18, 22, 5.0d);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(new TreeSet<>(filters.keySet()), fails, "Failed sites did not match expected set of failed sites.");
    }

    /** Utility method that takes a a VCF and a set of parameters and filters the VCF. */
    File testFiltering(final File vcf, final double minAb, final int minDp, final int minGq, final double maxFs) throws Exception {
        final File out = BaseTest.createTempFile("filterVcfTest.", ".vcf.gz");

        final String[] args = new String[]{
                "--CREATE_INDEX", "true",
                "--input", vcf.getAbsolutePath(),
                "--output", out.getAbsolutePath(),
                "--MIN_AB", Double.toString(minAb),
                "--MIN_DP", Integer.toString(minDp),
                "--MIN_GQ", Integer.toString(minGq),
                "--MAX_FS", Double.toString(maxFs)
        };

        runCommandLine(args);
        return out;
    }

    /** Consumes a VCF and returns a ListMap where each they keys are the IDs of filtered out sites and the values are the set of filters. */
    ListMap<String,String> slurpFilters(final File vcf) {
        final ListMap<String,String> map = new ListMap<>();
        try (final VCFFileReader in = new VCFFileReader(vcf, false)) {
            for (final VariantContext ctx : in) {
                if (ctx.isNotFiltered()) continue;
                for (final String filter : ctx.getFilters()) {
                    map.add(ctx.getID(), filter);
                }
            }
        }
        return map;
    }
}
