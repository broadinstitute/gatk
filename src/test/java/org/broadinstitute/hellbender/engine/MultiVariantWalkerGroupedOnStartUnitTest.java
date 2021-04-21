package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections.IteratorUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.mockito.Mockito;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.mockito.Mockito.when;

public class MultiVariantWalkerGroupedOnStartUnitTest extends GATKBaseTest {

    public static final Allele REF = Allele.create("A", true);
    public static final Allele ALT = Allele.create("C", false);

    private final class DummyExampleGroupingMultiVariantWalker extends MultiVariantWalkerGroupedOnStart {
        List<List<VariantContext>> seenVariants = new ArrayList<>();
        @Override
        public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
            //NOTE: MultiVariantDataSource stores the FeatureInput name as the source, which is machine-specific, so remove it:
            seenVariants.add(variantContexts.stream().map(vc -> new VariantContextBuilder(vc).source("Unknown").make()).collect(Collectors.toList()));
        }
    }

    @DataProvider
    public Object[][] vcfExampleFiles() {
        return new Object[][]{
                {getTestFile("gvcfExample1.vcf"), getTestFile("gvcfExample1.copy.vcf")},
                {getTestFile("gvcfExample2.vcf"), getTestFile("gvcfExample2.copy.vcf")}
        };
    }

    @Test (dataProvider = "vcfExampleFiles")
    @SuppressWarnings({"unchecked"})
    public void testSimpleCaseDuplicateFile(File file, File file2) {
        final DummyExampleGroupingMultiVariantWalker tool = new DummyExampleGroupingMultiVariantWalker();
        final String[] args = {"--variant", file.getAbsolutePath(), "--variant", file2.getAbsolutePath(), "-R", b37_reference_20_21};
        tool.instanceMain(args);

        //Check that we get the right number of apply calls and the right number of variants
        try(final FeatureDataSource<VariantContext> variantContextFeatureDataSource = new FeatureDataSource<>(file)) {
            List<VariantContext> inputList = IteratorUtils.toList(variantContextFeatureDataSource.iterator());

            Assert.assertEquals(tool.seenVariants.size(), inputList.size());

            // We expect each variant in a duplicate file to be applied once but duplicated across the files
            for (int i = 0; i < inputList.size(); i++) {
                Assert.assertEquals(tool.seenVariants.get(i).size(),2);
                Assert.assertEquals(tool.seenVariants.get(i).get(0).getStart(),inputList.get(i).getStart());
                Assert.assertEquals(tool.seenVariants.get(i).get(1).getStart(),inputList.get(i).getStart());
            }
        }
    }

    @Test
    public void testIgnoreVariantsThatStartOutsideInterval() {
        final DummyExampleGroupingMultiVariantWalker tool = new DummyExampleGroupingMultiVariantWalker();
        final String[] args = {"--variant", getTestFile("gvcfExample1.vcf").getAbsolutePath(), "--"+MultiVariantWalkerGroupedOnStart.IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL, "-L", "20:69500-69540", "-R", b37_reference_20_21};
        tool.instanceMain(args);

        Assert.assertEquals(tool.seenVariants.size(), 3);
        Assert.assertEquals(tool.seenVariants.get(0).get(0).getStart(), 69511); //Asserting that we skipped the first read which was partially in the window
        Assert.assertEquals(tool.seenVariants.get(2).get(0).getStart(), 69522); //Asserting that the last read which ends after the region is still present
    }

    @Test
    @SuppressWarnings({"unchecked"})
    public void testMergingTwoFilesCorrectBehavior() {
        final DummyExampleGroupingMultiVariantWalker tool = new DummyExampleGroupingMultiVariantWalker();
        final String[] args = {"--variant", getToolTestDataDir()+"gvcfExample1.vcf", "--variant", getToolTestDataDir()+"gvcfExample2.vcf", "-R", b37_reference_20_21};
        tool.instanceMain(args);

        //Check that we get the right number of apply calls and the right number of variants
        try(final FeatureDataSource<VariantContext> variantContextFeatureDataSource1 = new FeatureDataSource<>(getTestFile("gvcfExample1.vcf"));
            final FeatureDataSource<VariantContext> variantContextFeatureDataSource2 = new FeatureDataSource<>(getTestFile("gvcfExample2.vcf"))){
            List<String> inputList1 = ((List<VariantContext>)IteratorUtils.toList(variantContextFeatureDataSource1.iterator())).stream().map(VariantContext::toString).collect(Collectors.toList());
            List<String> inputList2 = ((List<VariantContext>)IteratorUtils.toList(variantContextFeatureDataSource2.iterator())).stream().map(VariantContext::toString).collect(Collectors.toList());

            Assert.assertEquals(tool.seenVariants.size(), 13);
            int total = 0;
            // We expect each variant in a duplicate file to be applied once but duplicated across the files
            for (int i = 0; i < tool.seenVariants.size(); i++) {
                total += tool.seenVariants.get(i).size();
                Assert.assertTrue(tool.seenVariants.get(i).size()==1 || tool.seenVariants.get(i).size()==2);
                Assert.assertTrue(inputList1.contains(tool.seenVariants.get(i).get(0).toString()) || inputList2.contains(tool.seenVariants.get(i).get(0).toString()));
            }
            Assert.assertEquals(total, inputList1.size()+inputList2.size());
        }
    }

    @Test
    @SuppressWarnings({"unchecked"})
    public void testSingleFileMultisites() {
        final DummyExampleGroupingMultiVariantWalker tool = new DummyExampleGroupingMultiVariantWalker();
        final String[] args = {"--variant", getToolTestDataDir()+"gvcfExample3.vcf", "-R", b37_reference_20_21};
        tool.instanceMain(args);

        //Check that we get the right number of apply calls and the right number of variants
        Assert.assertEquals(tool.seenVariants.size(), 1);
        // We expect all the variants piled up on one spot to be served at once
        Assert.assertEquals(tool.seenVariants.get(0).size(),7);
        List<String> variants = tool.seenVariants.get(0).stream().map(VariantContext::toString).collect(Collectors.toList());
        Assert.assertTrue(variants.stream().noneMatch(s -> variants.indexOf(s)!=variants.lastIndexOf(s)));
    }

}