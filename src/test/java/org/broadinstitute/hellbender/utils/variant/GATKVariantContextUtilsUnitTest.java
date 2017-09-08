package org.broadinstitute.hellbender.utils.variant;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public final class GATKVariantContextUtilsUnitTest extends BaseTest {

    Allele Aref, T, C, G, Cref, ATC, ATCATC;
    Allele ATCATCT;
    Allele ATref;
    Allele Anoref;
    Allele GT;

    @BeforeClass
    public void setup() throws IOException {
        // alleles
        Aref = Allele.create("A", true);
        Cref = Allele.create("C", true);
        T = Allele.create("T");
        C = Allele.create("C");
        G = Allele.create("G");
        ATC = Allele.create("ATC");
        ATCATC = Allele.create("ATCATC");
        ATCATCT = Allele.create("ATCATCT");
        ATref = Allele.create("AT",true);
        Anoref = Allele.create("A",false);
        GT = Allele.create("GT",false);
    }

    private Genotype makeG(String sample, Allele a1, Allele a2, double log10pError, int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).log10PError(log10pError).PL(pls).make();
    }


    private Genotype makeG(String sample, Allele a1, Allele a2, double log10pError) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).log10PError(log10pError).make();
    }

    private VariantContext makeVC(String source, List<Allele> alleles) {
        return makeVC(source, alleles, null, null);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Genotype... g1) {
        return makeVC(source, alleles, Arrays.asList(g1));
    }

    private VariantContext makeVC(String source, List<Allele> alleles, String filter) {
        return makeVC(source, alleles, filter.equals(".") ? null : new LinkedHashSet<String>(Collections.singletonList(filter)));
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Set<String> filters) {
        return makeVC(source, alleles, null, filters);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Collection<Genotype> genotypes) {
        return makeVC(source, alleles, genotypes, null);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Collection<Genotype> genotypes, Set<String> filters) {
        int start = 10;
        int stop = start + alleles.get(0).length() - 1; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(genotypes).filters(filters).make();
    }

    @Test
    public void testHomozygousAlleleList() throws Exception {
        final List<Allele> alleles = GATKVariantContextUtils.homozygousAlleleList(T, 2);
        Assert.assertEquals(alleles, Arrays.asList(T, T));
    }

    // --------------------------------------------------------------------------------
    //
    // Test allele merging
    //
    // --------------------------------------------------------------------------------

    private class MergeAllelesTest extends TestDataProvider {
        List<List<Allele>> inputs;
        List<Allele> expected;

        @SafeVarargs
        @SuppressWarnings("varargs")
        private MergeAllelesTest(List<Allele>... arg) {
            super(MergeAllelesTest.class);
            LinkedList<List<Allele>> all = new LinkedList<>(Arrays.asList(arg));
            expected = all.pollLast();
            inputs = all;
        }

        public String toString() {
            return String.format("MergeAllelesTest input=%s expected=%s", inputs, expected);
        }
    }
    @DataProvider(name = "mergeAlleles")
    public Object[][] mergeAllelesData() {
        // first, do no harm
        new MergeAllelesTest(Collections.singletonList(Aref), Collections.singletonList(Aref));

        new MergeAllelesTest(Collections.singletonList(Aref), Collections.singletonList(Aref), Collections.singletonList(Aref));

        new MergeAllelesTest(Collections.singletonList(Aref), Arrays.asList(Aref, T), Arrays.asList(Aref, T));

        new MergeAllelesTest(Arrays.asList(Aref, C), Arrays.asList(Aref, T), Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, T),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, T, C)); // in order of appearence

        new MergeAllelesTest(Arrays.asList(Aref, C, T),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, C, T), Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, T, C), Arrays.asList(Aref, T, C));

        new MergeAllelesTest(Arrays.asList(Aref, T, C),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, T, C)); // in order of appearence

        new MergeAllelesTest(Arrays.asList(Aref),
                Arrays.asList(Aref, ATC),
                Arrays.asList(Aref, ATC));

        new MergeAllelesTest(Arrays.asList(Aref),
                Arrays.asList(Aref, ATC, ATCATC),
                Arrays.asList(Aref, ATC, ATCATC));

        // alleles in the order we see them
        new MergeAllelesTest(Arrays.asList(Aref, ATCATC),
                Arrays.asList(Aref, ATC, ATCATC),
                Arrays.asList(Aref, ATCATC, ATC));

        // same
        new MergeAllelesTest(Arrays.asList(Aref, ATC),
                Arrays.asList(Aref, ATCATC),
                Arrays.asList(Aref, ATC, ATCATC));

        new MergeAllelesTest(Arrays.asList(ATref, ATC, Anoref, G),
                Arrays.asList(Aref, ATCATC, G),
                Arrays.asList(ATref, ATC, Anoref, G, ATCATCT, GT));

        return MergeAllelesTest.getTests(MergeAllelesTest.class);
    }

    @Test( dataProvider = "mergeAlleles")
    public void testMergeAlleles(MergeAllelesTest cfg) {
        final List<VariantContext> inputs = new ArrayList<>();

        int i = 0;
        for ( final List<Allele> alleles : cfg.inputs ) {
            final String name = "vcf" + ++i;
            inputs.add(makeVC(name, alleles));
        }

        final List<String> priority = vcs2priority(inputs);

        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                inputs, priority,
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, "set", false, false);

        Assert.assertEquals(merged.getAlleles().size(),cfg.expected.size());
        Assert.assertEquals(new LinkedHashSet<>(merged.getAlleles()), new LinkedHashSet<>(cfg.expected));   //HACK this is a hack to get around a bug in the htsjdk.  The method returns a list with an unspecified order.
    }

    // --------------------------------------------------------------------------------
    //
    // Test rsID merging
    //
    // --------------------------------------------------------------------------------

    private class SimpleMergeRSIDTest extends TestDataProvider {
        List<String> inputs;
        String expected;

        private SimpleMergeRSIDTest(String... arg) {
            super(SimpleMergeRSIDTest.class);
            LinkedList<String> allStrings = new LinkedList<>(Arrays.asList(arg));
            expected = allStrings.pollLast();
            inputs = allStrings;
        }

        public String toString() {
            return String.format("SimpleMergeRSIDTest vc=%s expected=%s", inputs, expected);
        }
    }

    @DataProvider(name = "simplemergersiddata")
    public Object[][] createSimpleMergeRSIDData() {
        new SimpleMergeRSIDTest(".", ".");
        new SimpleMergeRSIDTest(".", ".", ".");
        new SimpleMergeRSIDTest("rs1", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs1", "rs1");
        new SimpleMergeRSIDTest(".", "rs1", "rs1");
        new SimpleMergeRSIDTest("rs1", ".", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs1,rs2");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs1", "rs1,rs2"); // duplicates
        new SimpleMergeRSIDTest("rs2", "rs1", "rs2,rs1");
        new SimpleMergeRSIDTest("rs2", "rs1", ".", "rs2,rs1");
        new SimpleMergeRSIDTest("rs2", ".", "rs1", "rs2,rs1");
        new SimpleMergeRSIDTest("rs1", ".", ".", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs3", "rs1,rs2,rs3");

        return SimpleMergeRSIDTest.getTests(SimpleMergeRSIDTest.class);
    }

    @Test(dataProvider = "simplemergersiddata")
    public void testRSIDMerge(SimpleMergeRSIDTest cfg) {
        VariantContext snpVC1 = makeVC("snpvc1", Arrays.asList(Aref, T));
        final List<VariantContext> inputs = cfg.inputs.stream()
                .map(id -> new VariantContextBuilder(snpVC1).id(id).make())
                .collect(Collectors.toList());

        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                inputs, null,
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.UNSORTED, false, false, "set", false, false);
        Assert.assertEquals(merged.getID(), cfg.expected);
    }

    // --------------------------------------------------------------------------------
    //
    // Test filtered merging
    //
    // --------------------------------------------------------------------------------

    private class MergeFilteredTest extends TestDataProvider {
        List<VariantContext> inputs;
        VariantContext expected;
        String setExpected;
        GATKVariantContextUtils.FilteredRecordMergeType type;


        private MergeFilteredTest(String name, VariantContext input1, VariantContext input2, VariantContext expected, String setExpected) {
            this(name, input1, input2, expected, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, setExpected);
        }

        private MergeFilteredTest(String name, VariantContext input1, VariantContext input2, VariantContext expected, GATKVariantContextUtils.FilteredRecordMergeType type, String setExpected) {
            super(MergeFilteredTest.class, name);
            LinkedList<VariantContext> all = new LinkedList<>(Arrays.asList(input1, input2));
            this.expected = expected;
            this.type = type;
            inputs = all;
            this.setExpected = setExpected;
        }

        public String toString() {
            return String.format("%s input=%s expected=%s", super.toString(), inputs, expected);
        }
    }

    @DataProvider(name = "mergeFiltered")
    public Object[][] mergeFilteredData() {
        new MergeFilteredTest("AllPass",
                makeVC("1", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                GATKVariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("noFilters",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "."),
                GATKVariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("oneFiltered",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "."),
                String.format("1-%s2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("onePassOneFail",
                makeVC("1", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                String.format("1-%s2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("AllFiltered",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "FAIL"),
                GATKVariantContextUtils.MERGE_FILTER_IN_ALL);

        // test ALL vs. ANY
        new MergeFilteredTest("FailOneUnfiltered",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "."),
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                String.format("%s1-2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("OneFailAllUnfilteredArg",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "FAIL"),
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ALL_UNFILTERED,
                String.format("%s1-2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        // test excluding allele in filtered record
        new MergeFilteredTest("DontIncludeAlleleOfFilteredRecords",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "."),
                String.format("1-%s2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        // promotion of site from unfiltered to PASSES
        new MergeFilteredTest("UnfilteredPlusPassIsPass",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                GATKVariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("RefInAll",
                makeVC("1", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                GATKVariantContextUtils.MERGE_REF_IN_ALL);

        new MergeFilteredTest("RefInOne",
                makeVC("1", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                "2");

        return MergeFilteredTest.getTests(MergeFilteredTest.class);
    }

    @Test(dataProvider = "mergeFiltered")
    public void testMergeFiltered(MergeFilteredTest cfg) {
        final List<String> priority = vcs2priority(cfg.inputs);
        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                cfg.inputs, priority, cfg.type, GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, true, false, "set", false, false);

        // test alleles are equal
        Assert.assertEquals(merged.getAlleles(), cfg.expected.getAlleles());

        // test set field
        Assert.assertEquals(merged.getAttribute("set"), cfg.setExpected);

        // test filter field
        Assert.assertEquals(merged.getFilters(), cfg.expected.getFilters());
    }

    @Test
    public void testEqualSites() throws Exception {
        final VariantContext v1 = new VariantContextBuilder("foo", "1", 10, 10, Collections.singletonList(Aref)).make();
        final VariantContext v2 = new VariantContextBuilder("foo", "1", 11, 11, Collections.singletonList(Aref)).make();
        final VariantContext v1WithAlleles = new VariantContextBuilder("foo", "1", 11, 11, Arrays.asList(T, Aref)).make();
        Assert.assertTrue(GATKVariantContextUtils.equalSites(v1, v1));
        Assert.assertTrue(GATKVariantContextUtils.equalSites(v2, v2));

        Assert.assertFalse(GATKVariantContextUtils.equalSites(v1, v2));
        Assert.assertFalse(GATKVariantContextUtils.equalSites(v1, v1WithAlleles));
        Assert.assertFalse(GATKVariantContextUtils.equalSites(v2, v1));
        Assert.assertFalse(GATKVariantContextUtils.equalSites(v1WithAlleles, v1));
    }

    // --------------------------------------------------------------------------------
    //
    // Test genotype merging
    //
    // --------------------------------------------------------------------------------

    private class MergeGenotypesTest extends TestDataProvider {
        List<VariantContext> inputs;
        VariantContext expected;
        List<String> priority;

        private MergeGenotypesTest(String name, String priority, VariantContext... arg) {
            super(MergeGenotypesTest.class, name);
            LinkedList<VariantContext> all = new LinkedList<>(Arrays.asList(arg));
            this.expected = all.pollLast();
            inputs = all;
            this.priority = Arrays.asList(priority.split(","));
        }

        public String toString() {
            return String.format("%s input=%s expected=%s", super.toString(), inputs, expected);
        }
    }

    @DataProvider(name = "mergeGenotypes")
    public Object[][] mergeGenotypesData() {
        new MergeGenotypesTest("TakeGenotypeByPriority-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)));

        new MergeGenotypesTest("TakeGenotypeByPriority-1,2-nocall", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, -1)));

        new MergeGenotypesTest("TakeGenotypeByPriority-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2)));

        new MergeGenotypesTest("NonOverlappingGenotypes", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s2", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1), makeG("s2", Aref, T, -2)));

        new MergeGenotypesTest("PreserveNoCall", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s2", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, -1), makeG("s2", Aref, T, -2)));

        new MergeGenotypesTest("PerserveAlleles", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, C), makeG("s2", Aref, C, -2)),
                makeVC("3", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1), makeG("s2", Aref, C, -2)));

        new MergeGenotypesTest("TakeGenotypePartialOverlap-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2), makeG("s3", Aref, T, -3)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1), makeG("s3", Aref, T, -3)));

        new MergeGenotypesTest("TakeGenotypePartialOverlap-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2), makeG("s3", Aref, T, -3)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2), makeG("s3", Aref, T, -3)));

        //
        // merging genothpes with PLs
        //

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs", "1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1, 1, 2, 3)),
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1, 1, 2, 3)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles", "1",
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles-2", "1",
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles-2", "1",
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s2", Aref, C, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6), makeG("s2", Aref, C, -1, 1, 2, 3, 4, 5, 6)));

        new MergeGenotypesTest("TakeGenotypePartialOverlapWithPLs-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1,5,0,3)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2,4,0,2), makeG("s3", Aref, T, -3,3,0,2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2,4,0,2), makeG("s3", Aref, T, -3,3,0,2)));

        new MergeGenotypesTest("TakeGenotypePartialOverlapWithPLs-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref,ATC), makeG("s1", Aref, ATC, -1,5,0,3)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2,4,0,2), makeG("s3", Aref, T, -3,3,0,2)),
                // no likelihoods on result since type changes to mixed multiallelic
                makeVC("3", Arrays.asList(Aref, ATC, T), makeG("s1", Aref, ATC, -1), makeG("s3", Aref, T, -3)));

        new MergeGenotypesTest("MultipleSamplePLsDifferentOrder", "1,2",
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, C, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("2", Arrays.asList(Aref, T, C), makeG("s2", Aref, T, -2, 6, 5, 4, 3, 2, 1)),
                // no likelihoods on result since type changes to mixed multiallelic
                makeVC("3", Arrays.asList(Aref, C, T), makeG("s1", Aref, C, -1), makeG("s2", Aref, T, -2)));

        return MergeGenotypesTest.getTests(MergeGenotypesTest.class);
    }

    @Test(dataProvider = "mergeGenotypes")
    public void testMergeGenotypes(MergeGenotypesTest cfg) {
        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                cfg.inputs, cfg.priority, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, true, false, "set", false, false);

        // test alleles are equal
        Assert.assertEquals(merged.getAlleles(), cfg.expected.getAlleles());

        // test genotypes
        assertGenotypesAreMostlyEqual(merged.getGenotypes(), cfg.expected.getGenotypes());
    }

    // necessary to not overload equals for genotypes
    private void assertGenotypesAreMostlyEqual(GenotypesContext actual, GenotypesContext expected) {
        if (actual == expected) {
            return;
        }

        if (actual == null || expected == null) {
            Assert.fail("Maps not equal: expected: " + expected + " and actual: " + actual);
        }

        if (actual.size() != expected.size()) {
            Assert.fail("Maps do not have the same size:" + actual.size() + " != " + expected.size());
        }

        for (Genotype value : actual) {
            Genotype expectedValue = expected.get(value.getSampleName());

            Assert.assertEquals(value.getAlleles(), expectedValue.getAlleles(), "Alleles in Genotype aren't equal");
            Assert.assertEquals(value.getGQ(), expectedValue.getGQ(), "GQ values aren't equal");
            Assert.assertEquals(value.hasLikelihoods(), expectedValue.hasLikelihoods(), "Either both have likelihoods or both not");
            if ( value.hasLikelihoods() )
                Assert.assertEquals(value.getLikelihoods().getAsVector(), expectedValue.getLikelihoods().getAsVector(), "Genotype likelihoods aren't equal");
        }
    }

    @Test
    public void testMergeGenotypesUniquify() {
        final VariantContext vc1 = makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1));
        final VariantContext vc2 = makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2));

        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                Arrays.asList(vc1, vc2), null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.UNIQUIFY, false, false, "set", false, false);

        // test genotypes
        Assert.assertEquals(merged.getSampleNames(), new LinkedHashSet<>(Arrays.asList("s1.1", "s1.2")));
    }

// TODO: remove after testing
//    @Test(expectedExceptions = IllegalStateException.class)
//    public void testMergeGenotypesRequireUnique() {
//        final VariantContext vc1 = makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1));
//        final VariantContext vc2 = makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2));
//
//        final VariantContext merged = VariantContextUtils.simpleMerge(
//                Arrays.asList(vc1, vc2), null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
//                GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE, false, false, "set", false, false);
//    }

    // --------------------------------------------------------------------------------
    //
    // Misc. tests
    //
    // --------------------------------------------------------------------------------

    @Test
    public void testAnnotationSet() {
        for ( final boolean annotate : Arrays.asList(true, false)) {
            for ( final String set : Arrays.asList("set", "combine", "x")) {
                final List<String> priority = Arrays.asList("1", "2");
                VariantContext vc1 = makeVC("1", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS);
                VariantContext vc2 = makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS);

                final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                        Arrays.asList(vc1, vc2), priority, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                        GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, annotate, false, set, false, false);

                if ( annotate )
                    Assert.assertEquals(merged.getAttribute(set), GATKVariantContextUtils.MERGE_INTERSECTION);
                else
                    Assert.assertFalse(merged.hasAttribute(set));
            }
        }
    }

    private static List<String> vcs2priority(final Collection<VariantContext> vcs) {
        return vcs.stream()
                .map(VariantContext::getSource)
                .collect(Collectors.toList());
    }

    // --------------------------------------------------------------------------------
    //
    // basic allele clipping test
    //
    // --------------------------------------------------------------------------------

    private class ReverseClippingPositionTestProvider extends TestDataProvider {
        final String ref;
        final List<Allele> alleles = new ArrayList<>();
        final int expectedClip;

        private ReverseClippingPositionTestProvider(final int expectedClip, final String ref, final String... alleles) {
            super(ReverseClippingPositionTestProvider.class);
            this.ref = ref;
            for ( final String allele : alleles )
                this.alleles.add(Allele.create(allele));
            this.expectedClip = expectedClip;
        }

        @Override
        public String toString() {
            return String.format("ref=%s allele=%s reverse clip %d", ref, alleles, expectedClip);
        }
    }

    @DataProvider(name = "ReverseClippingPositionTestProvider")
    public Object[][] makeReverseClippingPositionTestProvider() {
        // pair clipping
        new ReverseClippingPositionTestProvider(0, "ATT", "CCG");
        new ReverseClippingPositionTestProvider(1, "ATT", "CCT");
        new ReverseClippingPositionTestProvider(2, "ATT", "CTT");
        new ReverseClippingPositionTestProvider(2, "ATT", "ATT");  // cannot completely clip allele

        // triplets
        new ReverseClippingPositionTestProvider(0, "ATT", "CTT", "CGG");
        new ReverseClippingPositionTestProvider(1, "ATT", "CTT", "CGT"); // the T can go
        new ReverseClippingPositionTestProvider(2, "ATT", "CTT", "CTT"); // both Ts can go

        return ReverseClippingPositionTestProvider.getTests(ReverseClippingPositionTestProvider.class);
    }

    @Test(dataProvider = "ReverseClippingPositionTestProvider")
    public void testReverseClippingPositionTestProvider(ReverseClippingPositionTestProvider cfg) {
        int result = GATKVariantContextUtils.computeReverseClipping(cfg.alleles, cfg.ref.getBytes());
        Assert.assertEquals(result, cfg.expectedClip);
    }


    // --------------------------------------------------------------------------------
    //
    // test splitting into bi-allelics
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "SplitBiallelics")
    public Object[][] makeSplitBiallelics() throws CloneNotSupportedException {
        List<Object[]> tests = new ArrayList<>();

        final VariantContextBuilder root = new VariantContextBuilder("x", "20", 10, 10, Arrays.asList(Aref, C));

        // biallelic -> biallelic
        tests.add(new Object[]{root.make(), Collections.singletonList(root.make())});

        // monos -> monos
        root.alleles(Collections.singletonList(Aref));
        tests.add(new Object[]{root.make(), Collections.singletonList(root.make())});

        root.alleles(Arrays.asList(Aref, C, T));
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(Aref, C)).make(),
                        root.alleles(Arrays.asList(Aref, T)).make())});

        root.alleles(Arrays.asList(Aref, C, T, G));
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(Aref, C)).make(),
                        root.alleles(Arrays.asList(Aref, T)).make(),
                        root.alleles(Arrays.asList(Aref, G)).make())});

        final Allele C      = Allele.create("C");
        final Allele CA      = Allele.create("CA");
        final Allele CAA     = Allele.create("CAA");
        final Allele CAAAA   = Allele.create("CAAAA");
        final Allele CAAAAA  = Allele.create("CAAAAA");
        final Allele Cref      = Allele.create("C", true);
        final Allele CAref     = Allele.create("CA", true);
        final Allele CAAref    = Allele.create("CAA", true);
        final Allele CAAAref   = Allele.create("CAAA", true);

        root.alleles(Arrays.asList(Cref, CA, CAA));
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(Cref, CA)).make(),
                        root.alleles(Arrays.asList(Cref, CAA)).make())});

        root.alleles(Arrays.asList(CAAref, C, CA)).stop(12);
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(CAAref, C)).make(),
                        root.alleles(Arrays.asList(CAref, C)).stop(11).make())});

        root.alleles(Arrays.asList(CAAAref, C, CA, CAA)).stop(13);
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(CAAAref, C)).make(),
                        root.alleles(Arrays.asList(CAAref, C)).stop(12).make(),
                        root.alleles(Arrays.asList(CAref, C)).stop(11).make())});

        root.alleles(Arrays.asList(CAAAref, CAAAAA, CAAAA, CAA, C)).stop(13);
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(Cref, CAA)).stop(10).make(),
                        root.alleles(Arrays.asList(Cref, CA)).stop(10).make(),
                        root.alleles(Arrays.asList(CAref, C)).stop(11).make(),
                        root.alleles(Arrays.asList(CAAAref, C)).stop(13).make())});

        final Allele threeCopies = Allele.create("GTTTTATTTTATTTTA", true);
        final Allele twoCopies = Allele.create("GTTTTATTTTA", true);
        final Allele zeroCopies = Allele.create("G", false);
        final Allele oneCopies = Allele.create("GTTTTA", false);
        tests.add(new Object[]{root.alleles(Arrays.asList(threeCopies, zeroCopies, oneCopies)).stop(25).make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(threeCopies, zeroCopies)).stop(25).make(),
                        root.alleles(Arrays.asList(twoCopies, zeroCopies)).stop(20).make())});

        return tests.toArray(new Object[][]{});
    }

    // --------------------------------------------------------------------------------
    //
    // Test repeats
    //
    // --------------------------------------------------------------------------------

    private class RepeatDetectorTest extends TestDataProvider {
        String ref;
        boolean isTrueRepeat;
        VariantContext vc;

        private RepeatDetectorTest(boolean isTrueRepeat, String ref, String refAlleleString, String ... altAlleleStrings) {
            super(RepeatDetectorTest.class);
            this.isTrueRepeat = isTrueRepeat;
            this.ref = ref;

            List<Allele> alleles = new LinkedList<>();
            final Allele refAllele = Allele.create(refAlleleString, true);
            alleles.add(refAllele);
            for ( final String altString: altAlleleStrings) {
                final Allele alt = Allele.create(altString, false);
                alleles.add(alt);
            }

            VariantContextBuilder builder = new VariantContextBuilder("test", "chr1", 1, refAllele.length(), alleles);
            this.vc = builder.make();
        }

        public String toString() {
            return String.format("%s refBases=%s trueRepeat=%b vc=%s", super.toString(), ref, isTrueRepeat, vc);
        }
    }

    @DataProvider(name = "RepeatDetectorTest")
    public Object[][] makeRepeatDetectorTest() {
        new RepeatDetectorTest(true,  "NAAC", "N", "NA");
        new RepeatDetectorTest(true,  "NAAC", "NA", "N");
        new RepeatDetectorTest(false, "NAAC", "NAA", "N");
        new RepeatDetectorTest(false, "NAAC", "N", "NC");
        new RepeatDetectorTest(false, "AAC", "A", "C");

        // running out of ref bases => false
        new RepeatDetectorTest(false, "NAAC", "N", "NCAGTA");

        // complex repeats
        new RepeatDetectorTest(true,  "NATATATC", "N", "NAT");
        new RepeatDetectorTest(true,  "NATATATC", "N", "NATA");
        new RepeatDetectorTest(true,  "NATATATC", "N", "NATAT");
        new RepeatDetectorTest(true,  "NATATATC", "NAT", "N");
        new RepeatDetectorTest(false, "NATATATC", "NATA", "N");
        new RepeatDetectorTest(false, "NATATATC", "NATAT", "N");

        // multi-allelic
        new RepeatDetectorTest(true,  "NATATATC", "N", "NAT", "NATAT");
        new RepeatDetectorTest(true,  "NATATATC", "N", "NAT", "NATA");
        new RepeatDetectorTest(true,  "NATATATC", "NAT", "N", "NATAT");
        new RepeatDetectorTest(true,  "NATATATC", "NAT", "N", "NATA"); // two As
        new RepeatDetectorTest(false, "NATATATC", "NAT", "N", "NATC"); // false
        new RepeatDetectorTest(false, "NATATATC", "NAT", "N", "NCC"); // false
        new RepeatDetectorTest(false, "NATATATC", "NAT", "NATAT", "NCC"); // false

        return RepeatDetectorTest.getTests(RepeatDetectorTest.class);
    }

    @Test( dataProvider = "RepeatDetectorTest")
    public void testRepeatDetectorTest(RepeatDetectorTest cfg) {

        // test alleles are equal
        Assert.assertEquals(GATKVariantContextUtils.isTandemRepeat(cfg.vc, cfg.ref.getBytes()), cfg.isTrueRepeat);
    }

    @Test
    public void testFindNumberOfRepetitions() throws Exception {
        /*
         *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
         *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
         *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
         */
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "GATAT".getBytes(), false), 2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "GATAT".getBytes(), true), 0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("A".getBytes(), "ATATG".getBytes(), true), 1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "ATATG".getBytes(), true), 2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("CCC".getBytes(), "CCCCCCCC".getBytes(), true),2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("CCC".getBytes(), "CCCCCCCC".getBytes(), false),2);

        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("ATG".getBytes(), "ATGATGATGATG".getBytes(), true),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("G".getBytes(), "ATGATGATGATG".getBytes(), true),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("T".getBytes(), "T".getBytes(), true),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "ATGATGATCATG".getBytes(), true),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("CCCCCCCC".getBytes(), "CCC".getBytes(), true),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "AT".getBytes(), true), 1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "".getBytes(), true), 0); //empty test string

        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("ATG".getBytes(), "ATGATGATGATG".getBytes(), false),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("G".getBytes(), "ATGATGATGATG".getBytes(), false),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("T".getBytes(), "T".getBytes(), false),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "ATGATGATCATG".getBytes(), false),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("CCCCCCCC".getBytes(), "CCC".getBytes(), false),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "AT".getBytes(), false), 1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "".getBytes(), false), 0); //empty test string
    }

    @Test
    public void testFindNumberOfRepetitionsFullArray() throws Exception {
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("XXXATG".getBytes(),3 ,3, "ATGATGATGATGYYY".getBytes(), 0, 12, true),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("GGGG".getBytes(), 0 ,1 , "GGGGATGATGATGATG".getBytes(), 4, 12, true),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("T".getBytes(), 0, 1, "TTTTT".getBytes(), 0, 1, true),1);

        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), 0, 2, "AT".getBytes(), 0, 0, true), 0); //empty test string
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), 0, 2, "AT".getBytes(), 1, 0, true), 0); //empty test string
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), 0, 2, "".getBytes(), 0, 0, true), 0); //empty test string

        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("XXXAT".getBytes(), 3, 2, "XXXGATAT".getBytes(), 4, 4, false), 2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), 0, 2, "GATAT".getBytes(), 0, 5, false), 2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("ATG".getBytes(), 0, 3, "ATGATGATGATG".getBytes(), 0, 12, false),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("ATG".getBytes(), 0, 3, "ATGATGATGATGATG".getBytes(), 3, 12, false),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("G".getBytes(), 0, 1, "ATGATGATGATG".getBytes(), 0, 12, false),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("G".getBytes(), 0, 1, "ATGATGATGATGATG".getBytes(), 0, 12, false),1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFindNumberOfRepetitionsNullArg1() throws Exception {
        GATKVariantContextUtils.findNumberOfRepetitions(null, "AT".getBytes(), false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFindNumberOfRepetitionsNullArg2() throws Exception {
        GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), null, false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFindNumberOfRepetitionsEmptyArg1() throws Exception {
        GATKVariantContextUtils.findNumberOfRepetitions("".getBytes(), "AT".getBytes(), false);
    }

    @Test
    public void testRepeatAllele() {
        Allele nullR = Aref;
        Allele nullA = Allele.create("A", false);
        Allele atc   = Allele.create("AATC", false);
        Allele atcatc   = Allele.create("AATCATC", false);
        Allele ccccR = Allele.create("ACCCC", true);
        Allele cc   = Allele.create("ACC", false);
        Allele cccccc   = Allele.create("ACCCCCC", false);
        Allele gagaR   = Allele.create("AGAGA", true);
        Allele gagagaga   = Allele.create("AGAGAGAGA", false);

        // - / ATC [ref] from 20-22
        String delLoc = "chr1";

        // - [ref] / ATC from 20-20
        String insLoc = "chr1";
        int insLocStart = 20;
        int insLocStop = 20;

        Pair<List<Integer>,byte[]> result;
        byte[] refBytes = "TATCATCATCGGA".getBytes();

        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("ATG".getBytes()),3);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("AAA".getBytes()),1);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("CACACAC".getBytes()),7);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("CACACA".getBytes()),2);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("CATGCATG".getBytes()),4);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("AATAATA".getBytes()),7);


        // A*,ATC, context = ATC ATC ATC : (ATC)3 -> (ATC)4
        VariantContext vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStop, Arrays.asList(nullR,atc)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],3);
        Assert.assertEquals(result.getLeft().toArray()[1],4);
        Assert.assertEquals(result.getRight().length,3);

        // ATC*,A,ATCATC
        vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStart+3, Arrays.asList(Allele.create("AATC", true),nullA,atcatc)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],3);
        Assert.assertEquals(result.getLeft().toArray()[1],2);
        Assert.assertEquals(result.getLeft().toArray()[2],4);
        Assert.assertEquals(result.getRight().length,3);

        // simple non-tandem deletion: CCCC*, -
        refBytes = "TCCCCCCCCATG".getBytes();
        vc = new VariantContextBuilder("foo", delLoc, 10, 14, Arrays.asList(ccccR,nullA)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],8);
        Assert.assertEquals(result.getLeft().toArray()[1],4);
        Assert.assertEquals(result.getRight().length,1);

        // CCCC*,CC,-,CCCCCC, context = CCC: (C)7 -> (C)5,(C)3,(C)9
        refBytes = "TCCCCCCCAGAGAGAG".getBytes();
        vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStart+4, Arrays.asList(ccccR,cc, nullA,cccccc)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],7);
        Assert.assertEquals(result.getLeft().toArray()[1],5);
        Assert.assertEquals(result.getLeft().toArray()[2],3);
        Assert.assertEquals(result.getLeft().toArray()[3],9);
        Assert.assertEquals(result.getRight().length,1);

        // GAGA*,-,GAGAGAGA
        refBytes = "TGAGAGAGAGATTT".getBytes();
        vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStart+4, Arrays.asList(gagaR, nullA,gagagaga)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],5);
        Assert.assertEquals(result.getLeft().toArray()[1],3);
        Assert.assertEquals(result.getLeft().toArray()[2],7);
        Assert.assertEquals(result.getRight().length,2);

    }

    // --------------------------------------------------------------------------------
    //
    // test forward clipping
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "ForwardClippingData")
    public Object[][] makeForwardClippingData() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{Collections.singletonList("A"), -1});
        tests.add(new Object[]{Collections.singletonList("<DEL>"), -1});
        tests.add(new Object[]{Arrays.asList("A", "C"), -1});
        tests.add(new Object[]{Arrays.asList("AC", "C"), -1});
        tests.add(new Object[]{Arrays.asList("A", "G"), -1});
        tests.add(new Object[]{Arrays.asList("A", "T"), -1});
        tests.add(new Object[]{Arrays.asList("GT", "CA"), -1});
        tests.add(new Object[]{Arrays.asList("GT", "CT"), -1});
        tests.add(new Object[]{Arrays.asList("ACC", "AC"), 0});
        tests.add(new Object[]{Arrays.asList("ACGC", "ACG"), 1});
        tests.add(new Object[]{Arrays.asList("ACGC", "ACG"), 1});
        tests.add(new Object[]{Arrays.asList("ACGC", "ACGA"), 2});
        tests.add(new Object[]{Arrays.asList("ACGC", "AGC"), 0});
        tests.add(new Object[]{Arrays.asList("A", "<DEL>"), -1});
        for ( int len = 0; len < 50; len++ )
            tests.add(new Object[]{Arrays.asList("A" + new String(Utils.dupBytes((byte)'C', len)), "C"), -1});

        tests.add(new Object[]{Arrays.asList("A", "T", "C"), -1});
        tests.add(new Object[]{Arrays.asList("AT", "AC", "AG"), 0});
        tests.add(new Object[]{Arrays.asList("AT", "AC", "A"), -1});
        tests.add(new Object[]{Arrays.asList("AT", "AC", "ACG"), 0});
        tests.add(new Object[]{Arrays.asList("AC", "AC", "ACG"), 0});
        tests.add(new Object[]{Arrays.asList("AC", "ACT", "ACG"), 0});
        tests.add(new Object[]{Arrays.asList("ACG", "ACGT", "ACGTA"), 1});
        tests.add(new Object[]{Arrays.asList("ACG", "ACGT", "ACGCA"), 1});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ForwardClippingData")
    public void testForwardClipping(final List<String> alleleStrings, final int expectedClip) {
        final List<Allele> alleles = alleleStrings.stream()
                .map(Allele::create)
                .collect(Collectors.toList());

        for ( final List<Allele> myAlleles : Utils.makePermutations(alleles, alleles.size(), false)) {
            final int actual = GATKVariantContextUtils.computeForwardClipping(myAlleles);
            Assert.assertEquals(actual, expectedClip);
        }
    }

    @DataProvider(name = "ClipAlleleTest")
    public Object[][] makeClipAlleleTest() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{Arrays.asList("ACC", "AC"), Arrays.asList("AC", "A"), 0});
        tests.add(new Object[]{Arrays.asList("ACGC", "ACG"), Arrays.asList("GC", "G"), 2});
        tests.add(new Object[]{Arrays.asList("ACGC", "ACGA"), Arrays.asList("C", "A"), 3});
        tests.add(new Object[]{Arrays.asList("ACGC", "AGC"), Arrays.asList("AC", "A"), 0});
        tests.add(new Object[]{Arrays.asList("AT", "AC", "AG"), Arrays.asList("T", "C", "G"), 1});
        tests.add(new Object[]{Arrays.asList("AT", "AC", "ACG"), Arrays.asList("T", "C", "CG"), 1});
        tests.add(new Object[]{Arrays.asList("AC", "ACT", "ACG"), Arrays.asList("C", "CT", "CG"), 1});
        tests.add(new Object[]{Arrays.asList("ACG", "ACGT", "ACGTA"), Arrays.asList("G", "GT", "GTA"), 2});
        tests.add(new Object[]{Arrays.asList("ACG", "ACGT", "ACGCA"), Arrays.asList("G", "GT", "GCA"), 2});

        // trims from left and right
        tests.add(new Object[]{Arrays.asList("ACGTT", "ACCTT"), Arrays.asList("G", "C"), 2});
        tests.add(new Object[]{Arrays.asList("ACGTT", "ACCCTT"), Arrays.asList("G", "CC"), 2});
        tests.add(new Object[]{Arrays.asList("ACGTT", "ACGCTT"), Arrays.asList("G", "GC"), 2});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ClipAlleleTest")
    public void testClipAlleles(final List<String> alleleStrings, final List<String> expected, final int numLeftClipped) {
        final int start = 10;
        final VariantContext unclipped = GATKVariantContextUtils.makeFromAlleles("test", "20", start, alleleStrings);
        final VariantContext clipped = GATKVariantContextUtils.trimAlleles(unclipped, true, true);

        Assert.assertEquals(clipped.getStart(), unclipped.getStart() + numLeftClipped);
        for ( int i = 0; i < unclipped.getAlleles().size(); i++ ) {
            final Allele trimmed = clipped.getAlleles().get(i);
            Assert.assertEquals(trimmed.getBaseString(), expected.get(i));
        }
    }

    // --------------------------------------------------------------------------------
    //
    // test primitive allele splitting
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "PrimitiveAlleleSplittingData")
    public Object[][] makePrimitiveAlleleSplittingData() {
        List<Object[]> tests = new ArrayList<>();

        // no split
        tests.add(new Object[]{"A", "C", 0, null});
        tests.add(new Object[]{"A", "AC", 0, null});
        tests.add(new Object[]{"AC", "A", 0, null});

        // one split
        tests.add(new Object[]{"ACA", "GCA", 1, Collections.singletonList(0)});
        tests.add(new Object[]{"ACA", "AGA", 1, Collections.singletonList(1)});
        tests.add(new Object[]{"ACA", "ACG", 1, Collections.singletonList(2)});

        // two splits
        tests.add(new Object[]{"ACA", "GGA", 2, Arrays.asList(0, 1)});
        tests.add(new Object[]{"ACA", "GCG", 2, Arrays.asList(0, 2)});
        tests.add(new Object[]{"ACA", "AGG", 2, Arrays.asList(1, 2)});

        // three splits
        tests.add(new Object[]{"ACA", "GGG", 3, Arrays.asList(0, 1, 2)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PrimitiveAlleleSplittingData")
    public void testPrimitiveAlleleSplitting(final String ref, final String alt, final int expectedSplit, final List<Integer> variantPositions) {

        final int start = 10;
        final VariantContext vc = GATKVariantContextUtils.makeFromAlleles("test", "20", start, Arrays.asList(ref, alt));

        final List<VariantContext> result = GATKVariantContextUtils.splitIntoPrimitiveAlleles(vc);

        if ( expectedSplit > 0 ) {
            Assert.assertEquals(result.size(), expectedSplit);
            for ( int i = 0; i < variantPositions.size(); i++ ) {
                Assert.assertEquals(result.get(i).getStart(), start + variantPositions.get(i));
            }
        } else {
            Assert.assertEquals(result.size(), 1);
            Assert.assertEquals(vc, result.get(0));
        }
    }

    // --------------------------------------------------------------------------------
    //
    // test allele remapping
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "AlleleRemappingData")
    public Object[][] makeAlleleRemappingData() {
        List<Object[]> tests = new ArrayList<>();

        final Allele originalBase1 = Allele.create((byte)'A');
        final Allele originalBase2 = Allele.create((byte)'T');

        for ( final byte base1 : BaseUtils.BASES ) {
            for ( final byte base2 : BaseUtils.BASES ) {
                for ( final int numGenotypes : Arrays.asList(0, 1, 2, 5) ) {
                    Map<Allele, Allele> map = new LinkedHashMap<>(2);
                    map.put(originalBase1, Allele.create(base1));
                    map.put(originalBase2, Allele.create(base2));

                    tests.add(new Object[]{map, numGenotypes});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AlleleRemappingData")
    public void testAlleleRemapping(final Map<Allele, Allele> alleleMap, final int numGenotypes) {

        final GATKVariantContextUtils.AlleleMapper alleleMapper = new GATKVariantContextUtils.AlleleMapper(alleleMap);

        @SuppressWarnings("unchecked")
        final GenotypesContext originalGC = createGenotypesContext(numGenotypes, new ArrayList<>(alleleMap.keySet()));

        final GenotypesContext remappedGC = GATKVariantContextUtils.updateGenotypesWithMappedAlleles(originalGC, alleleMapper);

        for ( int i = 0; i < numGenotypes; i++ ) {
            final Genotype originalG = originalGC.get(String.format("%d", i));
            final Genotype remappedG = remappedGC.get(String.format("%d", i));

            Assert.assertEquals(originalG.getAlleles().size(), remappedG.getAlleles().size());
            for ( int j = 0; j < originalG.getAlleles().size(); j++ )
                Assert.assertEquals(remappedG.getAllele(j), alleleMap.get(originalG.getAllele(j)));
        }
    }

    private static GenotypesContext createGenotypesContext(final int numGenotypes, final List<Allele> alleles) {
        Utils.resetRandomGenerator();
        final Random random = Utils.getRandomGenerator();

        final GenotypesContext gc = GenotypesContext.create();
        for ( int i = 0; i < numGenotypes; i++ ) {
            // choose alleles at random
            final List<Allele> myAlleles = new ArrayList<>();
            myAlleles.add(alleles.get(random.nextInt(2)));
            myAlleles.add(alleles.get(random.nextInt(2)));

            final Genotype g = new GenotypeBuilder(String.format("%d", i)).alleles(myAlleles).make();
            gc.add(g);
        }

        return gc;
    }

    // --------------------------------------------------------------------------------
    //
    // Test subsetDiploidAlleles
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "SubsetAllelesData")
    public Object[][] makesubsetAllelesData() {
        List<Object[]> tests = new ArrayList<>();

        final List<Allele> AA = Arrays.asList(Aref,Aref);
        final List<Allele> AC = Arrays.asList(Aref,C);
        final List<Allele> CC = Arrays.asList(C,C);
        final List<Allele> AG = Arrays.asList(Aref,G);
        final List<Allele> GG = Arrays.asList(G,G);
        final List<Allele> ACG = Arrays.asList(Aref,C,G);

        final VariantContext vcBase = new VariantContextBuilder("test", "20", 10, 10, AC).make();

        // haploid, one alt allele
        final double[] haploidRefPL = MathUtils.normalizeFromRealSpace(new double[]{0.9, 0.1});
        final double[] haploidAltPL = MathUtils.normalizeFromRealSpace(new double[]{0.1, 0.9});
        final double[] haploidUninformative = new double[]{0, 0};

        // diploid, one alt allele
        final double[] homRefPL = MathUtils.normalizeFromRealSpace(new double[]{0.9, 0.09, 0.01});
        final double[] hetPL = MathUtils.normalizeFromRealSpace(new double[]{0.09, 0.9, 0.01});
        final double[] homVarPL = MathUtils.normalizeFromRealSpace(new double[]{0.01, 0.09, 0.9});
        final double[] uninformative = new double[]{0, 0, 0};

        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(50).make();

        // the simple case where no selection occurs
        final Genotype aHaploidGT = new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).AD(new int[]{10,2}).PL(haploidRefPL).GQ(8).make();
        final Genotype cHaploidGT = new GenotypeBuilder(base).alleles(Collections.singletonList(C)).AD(new int[]{10,2}).PL(haploidAltPL).GQ(8).make();
        final Genotype aaGT = new GenotypeBuilder(base).alleles(AA).AD(new int[]{10,2}).PL(homRefPL).GQ(8).make();
        final Genotype acGT = new GenotypeBuilder(base).alleles(AC).AD(new int[]{10,2}).PL(hetPL).GQ(8).make();
        final Genotype ccGT = new GenotypeBuilder(base).alleles(CC).AD(new int[]{10,2}).PL(homVarPL).GQ(8).make();

        // haploid
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(aHaploidGT).make(), AC, Collections.singletonList(new GenotypeBuilder(aHaploidGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(cHaploidGT).make(), AC, Collections.singletonList(new GenotypeBuilder(cHaploidGT).make())});
        // diploid
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(aaGT).make(), AC, Collections.singletonList(new GenotypeBuilder(aaGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(acGT).make(), AC, Collections.singletonList(new GenotypeBuilder(acGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(ccGT).make(), AC, Collections.singletonList(new GenotypeBuilder(ccGT).make())});

        // uninformative test cases
        // diploid
        final Genotype uninformativeGT = new GenotypeBuilder(base).alleles(CC).PL(uninformative).GQ(0).make();
        final Genotype emptyGT = new GenotypeBuilder(base).alleles(GATKVariantContextUtils.noCallAlleles(2)).noPL().noGQ().make();
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(uninformativeGT).make(), AC, Collections.singletonList(emptyGT)});
        // haploid
        final Genotype haploidUninformativeGT = new GenotypeBuilder(base).alleles(Collections.singletonList(C)).PL(haploidUninformative).GQ(0).make();
        final Genotype haplpoidEmptyGT = new GenotypeBuilder(base).alleles(GATKVariantContextUtils.noCallAlleles(1)).noPL().noGQ().make();
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(haploidUninformativeGT).make(), AC, Collections.singletonList(haplpoidEmptyGT)});

        // subsetting from 3 to 2 alleles
        // diploid PL order is: AA, AC, CC, AG, CG, GG (00, 01, 11, 02, 12, 22)
        final double[] homRef3AllelesPL = new double[]{0, -10, -20, -30, -40, -50};
        final double[] hetRefC3AllelesPL = new double[]{-10, 0, -20, -30, -40, -50};
        final double[] homC3AllelesPL = new double[]{-20, -10, 0, -30, -40, -50};
        final double[] hetRefG3AllelesPL = new double[]{-20, -10, -30, 0, -40, -50};
        final double[] hetCG3AllelesPL = new double[]{-20, -10, -30, -40, 0, -50};
        final double[] homG3AllelesPL = new double[]{-20, -10, -30, -40, -50, 0};
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(homRef3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{0, -10, -20}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(hetRefC3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(AC).PL(new double[]{-10, 0, -20}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(homC3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(CC).PL(new double[]{-20, -10, 0}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(hetRefG3AllelesPL).make()).make(),
                AG,
                Collections.singletonList(new GenotypeBuilder(base).alleles(AG).PL(new double[]{-20, 0, -50}).GQ(200).make())});
        // wow, scary -- bad output but discussed with Eric and we think this is the only thing that can be done
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(hetCG3AllelesPL).make()).make(),
                AG,
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{0, -20, -30}).GQ(200).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(homG3AllelesPL).make()).make(),
                AG,
                Collections.singletonList(new GenotypeBuilder(base).alleles(GG).PL(new double[]{-20, -40, 0}).GQ(200).make())});

        // haploid PL order is: A, C, G (0, 1, 2)
        final double[] haploidRef3AllelesPL = new double[]{0, -10, -20};
        final double[] haploidAltC3AllelesPL = new double[]{-10, 0, -20};
        final double[] haploidAltG3AllelesPL = new double[]{-20, -10, 0};
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).PL(haploidRef3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).PL(new double[]{0, -10}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(C)).PL(haploidAltC3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(C)).PL(new double[]{-10, 0}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(G)).PL(haploidAltG3AllelesPL).make()).make(),
                AG,
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(G)).PL(new double[]{-20, 0}).GQ(200).make())});

        return tests.toArray(new Object[][]{});
    }

    /*@Test(dataProvider = "SubsetAllelesData")
    public void testSubsetAlleles(final VariantContext inputVC,
                                  final List<Allele> allelesToUse,
                                  final List<Genotype> expectedGenotypes) {
        // initialize cache of allele anyploid indices
        for (final Genotype genotype : inputVC.getGenotypes()) {
            GenotypeLikelihoods.initializeAnyploidPLIndexToAlleleIndices(inputVC.getNAlleles() - 1, genotype.getPloidy());
        }

        final GenotypesContext actual = AlleleSubsettingUtils.subsetAlleles(inputVC, allelesToUse, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN);

        Assert.assertEquals(actual.size(), expectedGenotypes.size());
        for ( final Genotype expected : expectedGenotypes ) {
            final Genotype actualGT = actual.get(expected.getSampleName());
            Assert.assertNotNull(actualGT);
            VariantContextTestUtils.assertGenotypesAreEqual(actualGT, expected);
        }
    }*/

    @DataProvider(name = "MakeGenotypeCallData")
    public Object[][] makeGenotypeCallData() {
        final List<Object[]> tests = new ArrayList<>();

        final List<Allele> AA = Arrays.asList(Aref,Aref);
        final List<Allele> AC = Arrays.asList(Aref,C);
        final List<Allele> CC = Arrays.asList(C,C);
        final List<Allele> AG = Arrays.asList(Aref,G);
        final List<Allele> CG = Arrays.asList(C,G);
        final List<Allele> GG = Arrays.asList(G,G);
        final List<Allele> AAA = Arrays.asList(Aref,Aref,Aref);
        final List<Allele> AAC = Arrays.asList(Aref,Aref,C);
        final List<Allele> ACC = Arrays.asList(Aref,C,C);
        final List<Allele> CCC = Arrays.asList(C,C,C);
        final List<Allele> AAG = Arrays.asList(Aref,Aref,G);
        final List<Allele> ACG = Arrays.asList(Aref,C,G);
        final List<Allele> CCG = Arrays.asList(C,C,G);
        final List<Allele> AGG = Arrays.asList(Aref,G,G);
        final List<Allele> CGG = Arrays.asList(C,G,G);
        final List<Allele> GGG = Arrays.asList(G,G,G);
        final List<List<Allele>> allDiploidSubsetAlleles = Arrays.asList(AC,AG,ACG);
        final List<List<Allele>> allTriploidSubsetAlleles = Arrays.asList(AAA,AAC,ACC,CCC,AAG,ACG,CCG,AGG,CGG,GGG);

        // for P=1, the index of the genotype a is a
        final double[] aRefPL = new double[]{0.9, 0.09, 0.01};
        final double[] cPL = new double[]{0.09, 0.9, 0.01};
        final double[] gPL = new double[]{0.01, 0.09, 0.9};
        final List<double[]> allHaploidPLs = Arrays.asList(aRefPL, cPL, gPL);
        final List<List<Allele>> allHaploidSubsetAlleles = Arrays.asList(Collections.singletonList(Aref), Collections.singletonList(G));

        // for P=2 and N=1, the ordering is 00,01,11
        final double[] homRefPL = new double[]{0.9, 0.09, 0.01};
        final double[] hetPL = new double[]{0.09, 0.9, 0.01};
        final double[] homVarPL = new double[]{0.01, 0.09, 0.9};
        final double[] uninformative = new double[]{0.33, 0.33, 0.33};
        final List<double[]> allDiploidPLs = Arrays.asList(homRefPL, hetPL, homVarPL, uninformative);

        // for P=3 and N=2, the ordering is 000, 001, 011, 111, 002, 012, 112, 022, 122, 222
        final double[] aaaPL = new double[]{0.9, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] aacPL = new double[]{0.01, 0.9, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] accPL = new double[]{0.01, 0.02, 0.9, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] cccPL = new double[]{0.01, 0.02, 0.03, 0.9, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] aagPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.9, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] acgPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.9, 0.06, 0.07, 0.08, 0.09};
        final double[] ccgPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.9, 0.07, 0.08, 0.09};
        final double[] aggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.9, 0.08, 0.09};
        final double[] cggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.9, 0.09};
        final double[] gggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.9};
        final double[] uninformativeTriploid = new double[]{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        final List<double[]> allTriploidPLs = Arrays.asList(homRefPL, hetPL, homVarPL, uninformativeTriploid);

        for ( final List<Allele> alleles : allHaploidSubsetAlleles ) {
            tests.add(new Object[]{1, GenotypeAssignmentMethod.SET_TO_NO_CALL, allHaploidPLs.get(0), Collections.singletonList(Aref), alleles, GATKVariantContextUtils.noCallAlleles(1)});
        }

        for ( final List<Allele> alleles : allDiploidSubsetAlleles ) {
            tests.add(new Object[]{2, GenotypeAssignmentMethod.SET_TO_NO_CALL, allDiploidPLs.get(0), AA, alleles, GATKVariantContextUtils.noCallAlleles(2)});
        }

        for ( final List<Allele> alleles : allTriploidSubsetAlleles ) {
            tests.add(new Object[]{3, GenotypeAssignmentMethod.SET_TO_NO_CALL, allTriploidPLs.get(0), AAA, alleles, GATKVariantContextUtils.noCallAlleles(3)});
        }

        final List<Allele> originalHaploidGT = Collections.singletonList(Aref);
        final List<Allele> haploidAllelesToUse = Arrays.asList(Aref, C, G );
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aRefPL, originalHaploidGT, haploidAllelesToUse, Collections.singletonList(Aref)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cPL, originalHaploidGT, haploidAllelesToUse, Collections.singletonList(C)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, gPL, originalHaploidGT, haploidAllelesToUse, Collections.singletonList(G)});

        for ( final List<Allele> originalGT : Arrays.asList(AA, AC, CC, AG, CG, GG) ) {
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, homRefPL, originalGT, AC, AA});
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, hetPL, originalGT, AC, AC});
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, homVarPL, originalGT, AC, CC});
        }

        for ( final List<Allele> originalGT : allTriploidSubsetAlleles) {
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aaaPL, originalGT, ACG, AAA});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aacPL, originalGT, ACG, AAC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, accPL, originalGT, ACG, ACC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cccPL, originalGT, ACG, CCC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aagPL, originalGT, ACG, AAG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, acgPL, originalGT, ACG, ACG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, ccgPL, originalGT, ACG, CCG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aggPL, originalGT, ACG, AGG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cggPL, originalGT, ACG, CCG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, gggPL, originalGT, ACG, GGG});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MakeGenotypeCallData")
    public void testMakeGenotypeCall(final int ploidy,
                                     final GenotypeAssignmentMethod mode,
                                     final double[] likelihoods,
                                     final List<Allele> originalGT,
                                     final List<Allele> allelesToUse,
                                     final List<Allele> expectedAlleles) {
        Utils.validateArg(originalGT.size() == ploidy, "original call must be consistent with ploidy");

        final GenotypeBuilder gb = new GenotypeBuilder("test");
        final double[] logLikelhoods = MathUtils.normalizeLog10(likelihoods);

        GATKVariantContextUtils.makeGenotypeCall(originalGT.size(), gb, mode, logLikelhoods, allelesToUse);

        final Genotype g = gb.make();
        Assert.assertEquals(new LinkedHashSet<>(g.getAlleles()), new LinkedHashSet<>(expectedAlleles));
    }

    @Test()
    public void testSubsetToRef() {
        final Map<Genotype, Genotype> tests = new LinkedHashMap<>();

        for ( final List<Allele> alleles : Arrays.asList(Collections.singletonList(Aref), Collections.singletonList(C), Arrays.asList(Aref, C), Arrays.asList(Aref, C, C) ) ) {
            for ( final String name : Arrays.asList("test1", "test2") ) {
                final GenotypeBuilder builder = new GenotypeBuilder(name, alleles);
                builder.DP(10);
                builder.GQ(30);
                builder.AD(alleles.size() == 1 ? new int[]{1} : (alleles.size() == 2 ? new int[]{1, 2} : new int[]{1, 2, 3}));
                builder.PL(alleles.size() == 1 ? new int[]{1} : (alleles.size() == 2 ? new int[]{1, 2} : new int[]{1, 2, 3}));
                builder.attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY,
                        alleles.size() == 1 ? new int[]{1, 2}  : (alleles.size() == 2 ? new int[]{1, 2, 3, 4} : new int[]{1, 2, 3, 4, 5, 6}));
                final List<Allele> refs = Collections.nCopies(alleles.size(), Aref);
                tests.put(builder.make(), builder.alleles(refs).noAD().noPL().attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, null).make());
            }
        }

        for ( final int n : Arrays.asList(1, 2, 3) ) {
            for ( final List<Genotype> genotypes : Utils.makePermutations(new ArrayList<>(tests.keySet()), n, false) ) {
                final VariantContext vc = new VariantContextBuilder("test", "20", 1, 1, Arrays.asList(Aref, C)).genotypes(genotypes).make();
                final GenotypesContext gc = GATKVariantContextUtils.subsetToRefOnly(vc, 2);

                Assert.assertEquals(gc.size(), genotypes.size());
                for ( int i = 0; i < genotypes.size(); i++ ) {
                    VariantContextTestUtils.assertGenotypesAreEqual(gc.get(i), tests.get(genotypes.get(i)));
                }
            }
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Test methods for merging reference confidence VCs
    //
    // --------------------------------------------------------------------------------


    @Test(dataProvider = "indexOfAlleleData")
    public void testIndexOfAllele(final Allele reference, final List<Allele> altAlleles, final List<Allele> otherAlleles) {
        final List<Allele> alleles = new ArrayList<>(altAlleles.size() + 1);
        alleles.add(reference);
        alleles.addAll(altAlleles);
        final VariantContext vc = makeVC("Source", alleles);

        for (int i = 0; i < alleles.size(); i++) {
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,alleles.get(i),true,true,true),i);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc, alleles.get(i), false, true, true),i);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,alleles.get(i),true,true,false),i);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,alleles.get(i),false,true,false),i);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,Allele.create(alleles.get(i),true),true,true,true),i);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,Allele.create(alleles.get(i),true),true,true,false),-1);
            if (i == 0) {
                Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,alleles.get(i),true,false,true),-1);
                Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,alleles.get(i),false,false,true),-1);
                Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,alleles.get(i),true,false,false),-1);
                Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,alleles.get(i),false,false,false),-1);
                Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,Allele.create(alleles.get(i).getBases(),true),false,true,true),i);
                Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,Allele.create(alleles.get(i).getBases(),false),false,true,true),-1);
            } else {
                Assert.assertEquals(GATKVariantContextUtils.indexOfAltAllele(vc,alleles.get(i),true),i - 1);
                Assert.assertEquals(GATKVariantContextUtils.indexOfAltAllele(vc,alleles.get(i),false), i - 1);
                Assert.assertEquals(GATKVariantContextUtils.indexOfAltAllele(vc,Allele.create(alleles.get(i),true),true),i-1);
                Assert.assertEquals(GATKVariantContextUtils.indexOfAltAllele(vc,Allele.create(alleles.get(i),true),false),-1);
            }
        }

        for (final Allele other : otherAlleles) {
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc, other, true, true, true), -1);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,other,false,true,true),-1);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,other,true,true,false),-1);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,other,false,true,false),-1);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,other,true,false,true),-1);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,other,false,false,true),-1);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc,other,true,false,false),-1);
            Assert.assertEquals(GATKVariantContextUtils.indexOfAllele(vc, other, false, false, false),-1);
        }
    }

    @DataProvider(name = "indexOfAlleleData")
    public Iterator<Object[]> indexOfAlleleData() {

        final Allele[] ALTERNATIVE_ALLELES = new Allele[] { T, C, G, ATC, ATCATC};

        final int lastMask = 0x1F;

        return new Iterator<Object[]>() {

            int nextMask = 0;

            @Override
            public boolean hasNext() {
                return nextMask <= lastMask;
            }

            @Override
            public Object[] next() {

                int mask = nextMask++;
                final List<Allele> includedAlleles = new ArrayList<>(5);
                final List<Allele> excludedAlleles = new ArrayList<>(5);
                for (int i = 0; i < ALTERNATIVE_ALLELES.length; i++) {
                    ((mask & 1) == 1 ? includedAlleles : excludedAlleles).add(ALTERNATIVE_ALLELES[i]);
                    mask >>= 1;
                }
                return new Object[] { Aref , includedAlleles, excludedAlleles};
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Test(dataProvider="overlapWithData")
    public void testOverlapsWith(final VariantContext vc, final GenomeLoc genomeLoc) {
        final boolean expected;

        if (genomeLoc.isUnmapped())
            expected = false;
        else if (vc.getStart() > genomeLoc.getStop())
            expected = false;
        else if (vc.getEnd() < genomeLoc.getStart())
            expected = false;
        else if (!vc.getContig().equals(genomeLoc.getContig()))
            expected = false;
        else
            expected = true;

        Assert.assertEquals(GATKVariantContextUtils.overlapsRegion(vc, genomeLoc), expected);
    }


    private final String[] OVERLAP_WITH_CHROMOSOMES =  { "1", "2" };
    private final int[] OVERLAP_WITH_EVENT_SIZES =  { -10, -1, 0, 1, 10 }; // 0 == SNP , -X xbp deletion, +X xbp insertion.
    private final int[] OVERLAP_WITH_EVENT_STARTS = { 1000, 1001,
            1005, 1010,
            1009, 1011,
            2000 };

    @DataProvider(name="overlapWithData")
    public Object[][] overlapWithData() {

        final int totalLocations = OVERLAP_WITH_CHROMOSOMES.length * OVERLAP_WITH_EVENT_SIZES.length * OVERLAP_WITH_EVENT_STARTS.length + 1;
        final int totalEvents = OVERLAP_WITH_CHROMOSOMES.length * OVERLAP_WITH_EVENT_SIZES.length * OVERLAP_WITH_EVENT_STARTS.length;
        final GenomeLoc[] locs = new GenomeLoc[totalLocations];
        final VariantContext[] events = new VariantContext[totalEvents];

        generateAllLocationsAndVariantContextCombinations(OVERLAP_WITH_CHROMOSOMES, OVERLAP_WITH_EVENT_SIZES,
                OVERLAP_WITH_EVENT_STARTS, locs, events);

        return generateAllParameterCombinationsForOverlapWithData(locs, events);
    }

    private Object[][] generateAllParameterCombinationsForOverlapWithData(GenomeLoc[] locs, VariantContext[] events) {
        final List<Object[]> result = new LinkedList<>();
        for (final GenomeLoc loc : locs)
            for (final VariantContext event : events)
                result.add(new Object[] { event , loc });

        return result.toArray(new Object[result.size()][]);
    }

    private void generateAllLocationsAndVariantContextCombinations(final String[] chrs, final int[] eventSizes,
                                                                   final int[] eventStarts, final GenomeLoc[] locs,
                                                                   final VariantContext[] events) {
        int nextIndex = 0;
        for (final String chr : chrs )
            for (final int size : eventSizes )
                for (final int starts : eventStarts ) {
                    locs[nextIndex] = hg19GenomeLocParser.createGenomeLoc(chr,starts,starts + Math.max(0,size));
                    events[nextIndex++] = new VariantContextBuilder().source("test").loc(chr,starts,starts + Math.max(0,size)).alleles(Arrays.asList(
                            Allele.create(randomBases(size <= 0 ? 1 : size + 1, true), true), Allele.create(randomBases(size < 0 ? -size + 1 : 1, false), false))).make();
                }

        locs[nextIndex++]  = GenomeLoc.UNMAPPED;
    }

    @Test(dataProvider = "totalPloidyData")
    public void testTotalPloidy(final int[] ploidies, final int defaultPloidy, final int expected) {
        final Genotype[] genotypes = new Genotype[ploidies.length];
        final List<Allele> vcAlleles = Arrays.asList(Aref,C);
        for (int i = 0; i < genotypes.length; i++)
            genotypes[i] = new GenotypeBuilder().alleles(GATKVariantContextUtils.noCallAlleles(ploidies[i])).make();
        final VariantContext vc = new VariantContextBuilder().chr("seq1").genotypes(genotypes).alleles(vcAlleles).make();
        Assert.assertEquals(GATKVariantContextUtils.totalPloidy(vc,defaultPloidy),expected," " + defaultPloidy + " " + Arrays.toString(ploidies));
    }

    @DataProvider(name="totalPloidyData")
    public Object[][] totalPloidyData() {
        final Random rdn = Utils.getRandomGenerator();
        final List<Object[]> resultList = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            final int sampleCount = rdn.nextInt(10);

            int expected = 0;
            final int defaultPloidy = rdn.nextInt(10) + 1;
            final int[] plodies = new int[sampleCount];
            for (int j = 0; j < sampleCount; j++) {
                plodies[j] = rdn.nextInt(10);
                expected += plodies[j] == 0 ? defaultPloidy : plodies[j];
            }
            resultList.add(new Object[] { plodies, defaultPloidy, expected });
        }
        return resultList.toArray(new Object[100][]);
    }

    private byte[] randomBases(final int length, final boolean reference) {
        final byte[] bases = new byte[length];
        bases[0] = (byte) (reference  ? 'A' : 'C');
        BaseUtils.fillWithRandomBases(bases, 1, bases.length);
        return bases;
    }

    @Test
    public void testCreateVariantContextWriterNoOptions() {
        final File tmpDir = createTempDir("createVCFTest");
        final File outFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + ".vcf");

        final VariantContextWriter vcw = GATKVariantContextUtils.createVCFWriter(outFile, null, false);
        vcw.close();

        final File outFileIndex = new File(outFile.getAbsolutePath() + ".idx");
        final File outFileMD5 = new File(outFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outFile.exists(), "No output file was not created");
        Assert.assertFalse(outFileIndex.exists(), "The createIndex argument was not honored");
        Assert.assertFalse(outFileMD5.exists(), "The createMD5 argument was not honored");
    }

    @DataProvider(name="createVCFWriterData")
    public Object[][] createVCFWriterData() {
        return new Object[][]{
                {".vcf", ".idx", true, true},
                {".vcf", ".idx", false, true},
                {".vcf", ".idx", true, false},

                {".bcf", ".idx", true, true},
                {".bcf", ".idx", false, true},
                {".bcf", ".idx", true, false},

                {".vcf.bgz", ".tbi", true, true},
                // AbstractFeatureReader fails to recognize this as block compressed unless it has an accompanying index
                // Is that correct behavior ?
                //".vcf.bgz", ".tbi", false, true},
                {".vcf.gz", ".tbi", false, true},
                {".vcf.bgz", ".tbi", true, false},

                // defaults to .vcf
                {".tmp", ".idx", false, true},
        };
    }
    @Test(dataProvider = "createVCFWriterData")
    public void testCreateVCFWriterWithOptions(
            final String outputExtension,
            final String indexExtension,
            final boolean createIndex,
            final boolean createMD5) throws IOException {

        final File tmpDir = createTempDir("createVCFTest");
        final File outputFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + outputExtension);

        Options options[] = createIndex ?
                new Options[] {Options.INDEX_ON_THE_FLY} :
                new Options[] {};
        try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                outputFile,
                makeSimpleSequenceDictionary(),
                createMD5,
                options)) {
            writeHeader(writer);
        }

        final File outFileIndex = new File(outputFile.getAbsolutePath() + indexExtension);
        final File outFileMD5 = new File(outputFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outputFile.exists(), "No output file was not created");
        Assert.assertEquals(outFileIndex.exists(), createIndex, "The createIndex argument was not honored");
        Assert.assertEquals(outFileMD5.exists(), createMD5, "The createMD5 argument was not honored");

        verifyFileType(outputFile, outputExtension);
    }

    // just make sure we can read the file with the corresponding codec
    private void verifyFileType(
            final File resultVCFFile,
            final String outputExtension) {
        final FeatureCodec<? extends Feature, ?> featureCodec = FeatureManager.getCodecForFile(resultVCFFile);

        if (outputExtension.equals(".vcf") ||
            outputExtension.equals(".vcf.bgz") ||
            outputExtension.equals(".vcf.gz") ||
            outputExtension.equals(".tmp"))
        {
            Assert.assertEquals(featureCodec.getClass(), VCFCodec.class,
                    "Wrong codec selected for file " + resultVCFFile.getAbsolutePath());
        }
        else if (outputExtension.equals(".bcf")) {
            Assert.assertEquals(featureCodec.getClass(), BCF2Codec.class,
                    "Wrong codec selected for file " + resultVCFFile.getAbsolutePath());
        }
        else {
            throw new IllegalArgumentException("Unknown file extension in createVCFWriter test validation");
        }
    }

    @DataProvider(name="createVCFWriterLenientData")
    public Object[][] createVCFWriterLenientData() {
        return new Object[][]{
                {".vcf", ".idx", true, true},
                {".vcf", ".idx", false, true},
                {".vcf", ".idx", true, false}
        };
    }

    @Test(dataProvider = "createVCFWriterLenientData")
    public void testCreateVCFWriterLenientTrue(
            final String outputExtension,
            final String indexExtension,
            final boolean createIndex,
            final boolean createMD5) throws IOException {

        final File tmpDir = createTempDir("createVCFTest");
        final File outputFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + outputExtension);

        Options options[] = createIndex ?
                new Options[] {Options.ALLOW_MISSING_FIELDS_IN_HEADER, Options.INDEX_ON_THE_FLY} :
                new Options[] {Options.ALLOW_MISSING_FIELDS_IN_HEADER};
        try (final VariantContextWriter vcw = GATKVariantContextUtils.createVCFWriter(
                outputFile,
                makeSimpleSequenceDictionary(),
                createMD5,
                options)) {
            writeHeader(vcw);
            writeBadVariant(vcw);  // verify leniency by writing a bogus attribute
        }

        final File outFileIndex = new File(outputFile.getAbsolutePath() + indexExtension);
        final File outFileMD5 = new File(outputFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outputFile.exists(), "No output file was not created");
        Assert.assertEquals(outFileIndex.exists(), createIndex, "The createIndex argument was not honored");
        Assert.assertEquals(outFileMD5.exists(), createMD5, "The createMD5 argument was not honored");
    }

    @Test(dataProvider = "createVCFWriterLenientData", expectedExceptions = IllegalStateException.class)
    public void testCreateVCFWriterLenientFalse(
            final String outputExtension,
            final String indexExtension, // unused
            final boolean createIndex,
            final boolean createMD5) throws IOException {

        final File tmpDir = createTempDir("createVCFTest");
        final File outputFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + outputExtension);

        Options options[] = createIndex ?
                new Options[] {Options.INDEX_ON_THE_FLY} :
                new Options[] {};
        try (final VariantContextWriter vcw = GATKVariantContextUtils.createVCFWriter(
                outputFile,
                makeSimpleSequenceDictionary(),
                createMD5,
                options)) {
            writeHeader(vcw);
            writeBadVariant(vcw); // write a bad attribute and throw...
        }
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCreateVariantContextWriterNoReference() {
        // should throw due to lack of reference
        final File outFile = createTempFile("testVCFWriter", ".vcf");
        final VariantContextWriter vcw =
                     GATKVariantContextUtils.createVCFWriter(
                             outFile,
                             null,
                             true,
                             Options.INDEX_ON_THE_FLY);
        vcw.close();
    }

    // This test verifies that we can write valid .gz/.tbi pair using a .gz input file with a large header.
    // Specifically, we want to make sure that index queries on the result return the first variants emitted into
    // the file, and that we don't encounter https://github.com/broadinstitute/gatk/issues/2821 and/or
    // https://github.com/broadinstitute/gatk/issues/2801.
    @Test
    public void testOnTheFlyTabixCreation() throws IOException {
        final File inputGZIPFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/8_mutect2_sorted.vcf.gz");
        final File outputGZIPFile = createTempFile("testOnTheFlyTabixCreation", ".vcf.gz");

        long recordCount = 0;
        try (final VariantContextWriter vcfWriter = GATKVariantContextUtils.createVCFWriter(
                 outputGZIPFile,
                 null,
                 false,
                 Options.INDEX_ON_THE_FLY);
            final FeatureReader<VariantContext> inputFileReader =
                     AbstractFeatureReader.getFeatureReader(inputGZIPFile.getAbsolutePath(), new VCFCodec(), false))
        {
            vcfWriter.writeHeader((VCFHeader)inputFileReader.getHeader());
            final Iterator<VariantContext> it = inputFileReader.iterator();
            while (it.hasNext()) {
                vcfWriter.add(it.next());
                recordCount++;
            }
        }

        // make sure we got a tabix index
        final File tabixIndexFile = new File(outputGZIPFile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
        Assert.assertTrue(tabixIndexFile.exists());
        Assert.assertTrue(tabixIndexFile.length() > 0);

        // verify the index via query
        final SimpleInterval queryInterval = new SimpleInterval("chr6:33414233-118314029");
        try (final FeatureReader<VariantContext> outputFileReader = AbstractFeatureReader.getFeatureReader(
                outputGZIPFile.getAbsolutePath(),
                new VCFCodec(),
                true)) // require index
        {
            final long actualCount = outputFileReader.query(
                    queryInterval.getContig(), queryInterval.getStart(), queryInterval.getEnd()).stream().count();
            Assert.assertEquals(actualCount, recordCount);
        }
    }

    private void writeHeader(final VariantContextWriter writer) {
        final Set<VCFHeaderLine> metaData = new HashSet<>();
        metaData.add(new VCFHeaderLine(
                VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));
        final VCFHeader vcfHeader = new VCFHeader(metaData, Collections.emptyList());
        vcfHeader.setSequenceDictionary(makeSimpleSequenceDictionary());
        writer.writeHeader(vcfHeader);
    }

    private void writeBadVariant(final VariantContextWriter writer) {
        //write a variant with a (bad) attribute that doesn't appear in the header to the output
        final VariantContextBuilder vcBuilder = new VariantContextBuilder("","chr1", 1, 1, Arrays.asList(Aref));
        vcBuilder.attribute("fake", new Object());
        final VariantContext vc = vcBuilder.make();
        writer.add(vc);
    }

    private SAMSequenceDictionary makeSimpleSequenceDictionary() {
        final SAMSequenceDictionary seqDictionary = new SAMSequenceDictionary();
        seqDictionary.addSequence(new SAMSequenceRecord("chr1", 10));
        return seqDictionary;
    }

    @DataProvider(name = "setFilteredGenotypeToNocallTest")
    public Object[][] makeSetFilteredGenotypeToNocallTest() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new ArrayList<>(Arrays.asList(Aref, C)), true, Arrays.asList(0), 0, Arrays.asList(0.5)});
        tests.add(new Object[]{new ArrayList<>(Arrays.asList(Aref, C)), false, Arrays.asList(1), 2, Arrays.asList(0.5)});
        tests.add(new Object[]{new ArrayList<>(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)), true, Arrays.asList(1), 2, Arrays.asList(0.5)});

        return tests.toArray(new Object[][]{});
    }

    private List<String> getGenotypeFilters(final VariantContext vc, final Genotype g) {
        final List<String> filters = new ArrayList<>();
        if (g.isFiltered()) {
            filters.add(g.getFilters());
        }

        return filters;
    }

    @Test(dataProvider="setFilteredGenotypeToNocallTest")
    public void testSetFilteredGenotypeToNocall(final List<Allele> genotypeAlleles, final boolean setFilteredGenotypesToNocall,
                                                final List<Integer> calledAltAlleles, final int calledAlleles, final List<Double> alleleFrequency) {
        final List<Allele> alleles = new ArrayList<>(Arrays.asList(Aref, C));
        final Genotype genotype = new GenotypeBuilder("NA12878").alleles(genotypeAlleles).filter("FILTER").make();
        final VariantContextBuilder builder = new VariantContextBuilder("test", "chr1", 1, Aref.length(), alleles).genotypes(genotype).
                attribute(VCFConstants.ALLELE_COUNT_KEY, new int[]{1} ).
                attribute(VCFConstants.ALLELE_NUMBER_KEY, 2).
                attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[]{0.5});
        GATKVariantContextUtils.setFilteredGenotypeToNocall(builder, builder.make(), setFilteredGenotypesToNocall, this::getGenotypeFilters);
        VariantContext vc = builder.make();
        if ( vc.hasAttribute((VCFConstants.ALLELE_COUNT_KEY)) ) {
            final int[] array = (int[]) vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY);
            Assert.assertEquals(vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY), calledAltAlleles.toArray());
        }
        if ( vc.hasAttribute((VCFConstants.ALLELE_NUMBER_KEY)) ) {
            Assert.assertEquals(vc.getAttribute(VCFConstants.ALLELE_NUMBER_KEY), calledAlleles);
        }
        if ( vc.hasAttribute((VCFConstants.ALLELE_FREQUENCY_KEY)) ) {
            Assert.assertEquals(vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY), alleleFrequency.toArray());
        }
    }

    @DataProvider(name="gqFromPLsData")
    public Object[][] gqFromPLsData() {
        return new Object[][]{
                {new int[]{0, 15}, 15},
                {new int[]{15, 0}, 15},
                {new int[]{0, 10, 20}, 10},
                {new int[]{20, 10, 0}, 10},
                {new int[]{0, 10, 20, 30, 40}, 10},
                {new int[]{30, 40, 20, 10, 0}, 10},
                {new int[]{-10, 20, 35}, 30},
                {new int[]{35, 40, -10, 15, 20}, 25},
                {new int[]{0, 10, 20, 30, 40, 50, 5}, 5},
                {new int[]{15, 15, 0, 5}, 5},
                {new int[]{15, 15, 0, 25}, 15},
                {new int[]{0, 15, 0, 25}, 0}
        };
    }

    @Test(dataProvider = "gqFromPLsData")
    public void testCalculateGQFromPLs(final int[] plValues, final int expectedGQ) {
        Assert.assertEquals(GATKVariantContextUtils.calculateGQFromPLs(plValues), expectedGQ);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCalculateGQFromShortPLArray() {
        final int[] plValues = new int[]{0};
        GATKVariantContextUtils.calculateGQFromPLs(plValues);
    }

    @Test
    public void testCreateFilterListWithAppend() {
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Cref);
        alleles.add(Allele.create("A", false));
        final VariantContextBuilder vcBuilder = new VariantContextBuilder("","chr1", 1, 1, alleles);
        vcBuilder.filters("F1", "F2", "Y");
        final VariantContext vc = vcBuilder.make();

        final String testFilterString = "TEST_FILTER";
        final List<String> filterResult = GATKVariantContextUtils.createFilterListWithAppend(vc, testFilterString);
        Assert.assertEquals(filterResult.size(), 4);
        Assert.assertTrue(vc.getFilters().stream().allMatch(f -> filterResult.contains(f)));
        Assert.assertTrue(filterResult.contains(testFilterString));

        // Check that it is in the proper position
        Assert.assertEquals(filterResult.get(2), testFilterString);
    }

    @Test
    public void testCreateFilterListWithAppendToEmpty() {
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Cref);
        alleles.add(Allele.create("A", false));
        final VariantContextBuilder vcBuilder = new VariantContextBuilder("","chr1", 1, 1, alleles);
        final VariantContext vc = vcBuilder.make();

        final String testFilterString = "TEST_FILTER";
        final List<String> filterResult = GATKVariantContextUtils.createFilterListWithAppend(vc, testFilterString);
        Assert.assertEquals(filterResult.size(), 1);
        Assert.assertTrue(filterResult.contains(testFilterString));
    }
}
