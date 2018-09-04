package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A class to aid implementation of end-end integration tests over artificial data intended to test allele specific annotation implementations,
 * As of 8/29/17, test output was matched against GATK3 implementations of combineGVCFs().
 *
 * Future tool authors who seek to extend this class need only implement getAnnotationToUse(), getRawKey(), and getKey()
 * to invoke a series of tests mimicking the behavior of GATK3 combineGVCFs and genotypeGVCFs. If the annotation is not present
 * in the source files test will accept at that site. If one seeks to add new annotations to these once can edit the listed
 * sites in the source files NA12878.AS.chr20snippet.g.vcf, NA12892.AS.chr20snippet.g.vcf, CombineGVCFs.output.vcf, and
 * GenotypeGVCFs.output.vcf to add new annotations.
 */
public abstract class ReducibleAnnotationBaseTest extends GATKBaseTest {

    @Override
    public String getToolTestDataDir() {
        return "src/test/resources/" + this.getClass().getPackage().getName().replace(".", "/") + "/";
    }

    @DataProvider
    public Object[][] interestingSitesCombineResults() {
        List<Object[]> tests = new ArrayList<>();

        try (
                FeatureDataSource<VariantContext> vcfA = new FeatureDataSource<>(getTestFile("NA12878.AS.chr20snippet.g.vcf"));
                FeatureDataSource<VariantContext> vcfB = new FeatureDataSource<>(getTestFile("NA12892.AS.chr20snippet.g.vcf"));
                FeatureDataSource<VariantContext> combineVCFOutput = new FeatureDataSource<>(getTestFile("CombineGVCFs.output.vcf"));
                FeatureDataSource<VariantContext> genotypeGVCFOutput = new FeatureDataSource<>(getTestFile("GenotypeGVCFs.output.vcf"));
        ) {
            // these are hand picked sites from the allele specific unit tests for combinegvcfs that triggered a combine in GATK3
            Integer[] interestingLocs = {10087820, 10433312, 10433322, 10433324, 10433326, 10433328, 10433382, 10433391, 10433468, 10433560, 10433575, 10433594, 10433955,
                    10434177, 10434384, 10435067, 10434258, 10436227, 10684106};
            List<SimpleInterval> intervals = Arrays.stream(interestingLocs).map(m -> new SimpleInterval("20", m, m)).collect(Collectors.toList());
            for (SimpleInterval loc : intervals) {
                VariantContext a = vcfA.query(loc).next();
                VariantContext b = vcfB.query(loc).next();
                VariantContext result = combineVCFOutput.query(loc).next();
                Iterator<VariantContext> query = genotypeGVCFOutput.query(loc);
                VariantContext genotyped = query.hasNext() ? query.next() : null;
                tests.add(new Object[]{Arrays.asList(a, b), result, genotyped});
            }

        }
        return tests.toArray(new Object[][]{});
    }

    /*
     * Methods that must be overridden in order for the automatic GATK3 combineGVCFs tests to be run.
     */
    /**
     * This method needs to return a list of the discoverable class names for any annotations that are to be tested in this framework.
     */
    protected abstract List<Annotation> getAnnotationsToUse();

    /**
     * This should return the raw key produced by the tested reducible annotation so the tool can assert similarity to GATK3 combineGVCFs output
     */
    protected abstract String getRawKey();

    /**
     * This should return the final key produced by the tested reducible annotation so the tool can assert similarity to GATK3 genotypeGVCFs output
     */
    protected abstract String getKey();


    // NOTE: this code is mimicking the behavior of GATK3 combineGVCFS insofar as it is important for the annotations
    @Test(dataProvider = "interestingSitesCombineResults")
    public void testCombineAnnotationGATK3Concordance(List<VariantContext> VCs, VariantContext result, VariantContext genotyped) throws Exception {
        VariantAnnotatorEngine annotatorEngine = new VariantAnnotatorEngine(getAnnotationsToUse(), null, Collections.emptyList(), false);
        ReferenceConfidenceVariantContextMerger merger = new ReferenceConfidenceVariantContextMerger(annotatorEngine, new VCFHeader());
        VariantContext merged = merger.merge(VCs, new SimpleInterval(result.getContig(), result.getStart(), result.getStart()), result.getReference().getBases()[0], false, false);
        Assert.assertTrue(VariantContextTestUtils.alleleSpecificAnnotationEquals(merged, result, getRawKey()));
    }

    // NOTE: this code is mimicking the behavior of GATK3 GenotypeGVCFs
    @Test(dataProvider = "interestingSitesCombineResults")
    public void testFinalizeAnnotationGATK3Concordance(List<VariantContext> VCs, VariantContext result, VariantContext genotyped) throws Exception {
        if (result == null || genotyped == null) {
            return;
        }

        VariantAnnotatorEngine annotatorEngine = new VariantAnnotatorEngine(getAnnotationsToUse(), null, Collections.emptyList(), false);
        final UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = new GenotypeCalculationArgumentCollection();
        GeneralPloidyFailOverAFCalculatorProvider calculatorProvider = new GeneralPloidyFailOverAFCalculatorProvider(uac.genotypeArgs);

        GenotypingEngine<?> genotypingEngine = new MinimalGenotypingEngine(uac, new IndexedSampleList(result.getSampleNamesOrderedByName()), calculatorProvider);
        genotypingEngine.setAnnotationEngine(annotatorEngine);
        GenotypeLikelihoodsCalculationModel model = result.getType() == VariantContext.Type.INDEL
                ? GenotypeLikelihoodsCalculationModel.INDEL
                : GenotypeLikelihoodsCalculationModel.SNP;
        VariantContext withGenotypes = genotypingEngine.calculateGenotypes(result, model, null);
        withGenotypes = new VariantContextBuilder(withGenotypes).attributes(result.getAttributes()).make();
        VariantContext finalized = annotatorEngine.finalizeAnnotations(withGenotypes, result);
        finalized =  GATKVariantContextUtils.reverseTrimAlleles(finalized);
        Assert.assertTrue(VariantContextTestUtils.alleleSpecificAnnotationEquals(finalized, genotyped, getKey()));
    }
}
