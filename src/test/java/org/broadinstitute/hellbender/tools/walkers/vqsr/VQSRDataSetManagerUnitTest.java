package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.*;
import joptsimple.internal.Strings;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.apache.commons.io.output.NullOutputStream.NULL_OUTPUT_STREAM;

public class VQSRDataSetManagerUnitTest extends BaseTest {

    private final String TRUE_SNPS = getTestDataDir() + "/"+ "vqsr.trueSnps.vcf";
    private final String TRUE_INDELS = getTestDataDir() + "/"+ "vqsr.trueIndels.vcf";

    @Test
    public void testEmpty(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        Assert.assertFalse(vtsm.hasTrainingSet(), "training set");
        Assert.assertFalse(vtsm.hasTruthSet(), "truth set");
    }


    private static class TestClass{
        @Argument(fullName = "resource", shortName = "resource")
        List<FeatureInput<VariantContext>> resource;
    }

    private static void makeAndAddTrainingSets(VQSRTrainingSetManager vtsm, String... args) {
        TestClass obj= new TestClass();
        final String argList = Strings.join(Arrays.asList(args), " ");
        new CommandLineParser(obj).parseArguments(new PrintStream(NULL_OUTPUT_STREAM), argList.split(" ", -1));
        vtsm.addDataSets(obj.resource);
    }

    @Test
    public void testOneFeatureInput(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String args = "--resource resource,known=true,prior=10.0:myFile";
        makeAndAddTrainingSets(vtsm, args);
        Assert.assertFalse(vtsm.hasTrainingSet(), "training set");
        Assert.assertFalse(vtsm.hasTruthSet(), "truth set");
    }

    @Test
    public void testFeatureInputWithTruthSetTrainingSet(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String[] args= {"--resource resource,known=true,prior=10.0:myFile",
                              "--resource resource,truth=true,training=true,prior=15.0:myFile3"};

        makeAndAddTrainingSets(vtsm, args);

        Assert.assertTrue(vtsm.hasTrainingSet(), "training set");
        Assert.assertTrue(vtsm.hasTruthSet(), "truth set");
    }

    @Test
    public void testFeatureInputWithTruthSetTrainingSet_twoFiles(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String[] args= {"--resource resource,known=true,prior=10.0:myFile",
                              "--resource resource,truth=true,prior=15.0:myFile2",
                              "--resource resource,training=true,prior=15.0:myFile3"};
        makeAndAddTrainingSets(vtsm, args);

        Assert.assertTrue(vtsm.hasTrainingSet(), "training set");
        Assert.assertTrue(vtsm.hasTruthSet(), "truth set");
    }

    @Test
    public void testFeatureInputWithOnlyTrainingSet(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String[] args= {"--resource resource,known=true,prior=10.0:myFile",
                              "--resource resource,training=true,prior=15.0:myFile2"};
        makeAndAddTrainingSets(vtsm, args);

        Assert.assertTrue(vtsm.hasTrainingSet(), "training set");
        Assert.assertFalse(vtsm.hasTruthSet(), "truth set");
    }

    @Test
    public void testFeatureInputWithOnlyTruthSet(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String[] args= {"--resource resource,known=true,prior=10.0:myFile",
                              "--resource resource,truth=true,prior=15.0:myFile2"};
        makeAndAddTrainingSets(vtsm, args);
        Assert.assertFalse(vtsm.hasTrainingSet(), "training set");
        Assert.assertTrue(vtsm.hasTruthSet(), "truth set");
    }

    @Test
    public void testFeatureInputConsensusSet(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String[] args= {"--resource resource,consensus=true,prior=10.0:myFile",
                "--resource resource,truth=true,prior=15.0:myFile2"};
        makeAndAddTrainingSets(vtsm, args);
        Assert.assertFalse(vtsm.hasTrainingSet(), "training set");
        Assert.assertTrue(vtsm.hasTruthSet(), "truth set");
    }

    @Test
    public void testFeatureInputBadSet(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String[] args= {"--resource resource,bad=true,prior=10.0:myFile",
                "--resource resource,truth=true,prior=15.0:myFile2"};
        makeAndAddTrainingSets(vtsm, args);
        Assert.assertFalse(vtsm.hasTrainingSet(), "training set");
        Assert.assertTrue(vtsm.hasTruthSet(), "truth set");
    }

    //----------------------------------------------------------
    private static class TestCLP extends CommandLineProgram{
        @Argument(fullName = "resource", shortName = "resource")
        private List<FeatureInput<VariantContext>> resource= new ArrayList<>();

        @Override
        protected Object doWork() {
            return null;
        }
    }

    private static void makeAndAddTrainingSetsForCLP(TestCLP clp, VQSRTrainingSetManager vtsm, String... args) {
        final String argList = Strings.join(Arrays.asList(args), " ");
        new CommandLineParser(clp).parseArguments(new PrintStream(NULL_OUTPUT_STREAM), argList.split(" ", -1));
        vtsm.addDataSets(clp.resource);
    }

    @Test
    public void testSNPinTruthData(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String args = "--resource indels,training=true,truth=true,prior=15.0:" + TRUE_SNPS;
        TestCLP clp = new TestCLP();
        makeAndAddTrainingSetsForCLP(clp, vtsm, args);
        FeatureManager fm = new FeatureManager(clp);

        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(50).make();
        final String source= "fred";

        int start = 55550;//1	55550	rs2949421	A	T

        final Allele A = Allele.create("A", true);
        final Allele T = Allele.create("T");
        final List<Allele> AT = Arrays.asList(A, T);
        final double[] hetPL = MathUtils.normalizeFromRealSpace(new double[]{0.09, 0.9, 0.01});
        final Genotype atGT = new GenotypeBuilder(base).alleles(AT).AD(new int[]{10,2}).PL(hetPL).GQ(8).make();

        final List<Allele> alleles = AT;
        int stop = start + alleles.get(0).length() - 1; // alleles.contains(ATC) ? start + 3 : start;
        final Genotype genotypes = atGT;
        final String filters = null;
        final VariantContext vcSNP = new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(genotypes).filters(filters).make();

        final FeatureContext fc = new FeatureContext(fm, new SimpleInterval("1", start,	stop));
        VariantDatum vd = new VariantDatum();
        vtsm.annotateDatum(fc, vcSNP, vd, false);
        Assert.assertEquals(15.0, vd.prior, 0.0, "prior");  //vd matches -> prior gets adjusted
        Assert.assertFalse(vd.isKnown, "isKnown");
        Assert.assertTrue(vd.atTruthSite, "atTruthSite");
        Assert.assertTrue(vd.atTrainingSite, "atTrainingSite");
        Assert.assertEquals(0, vd.consensusCount, "consensusCount");

        final Allele ACT = Allele.create("ACT");
        final List<Allele> A_AT = Arrays.asList(A, ACT);
        final Genotype genotypesIndel = new GenotypeBuilder(base).alleles(A_AT).AD(new int[]{10,2}).PL(hetPL).GQ(8).make();;
        final List<Allele> allelesIndel = A_AT;
        final VariantContext vcIndel = new VariantContextBuilder(source, "1", start, stop, allelesIndel).genotypes(genotypesIndel).filters(filters).make();
        VariantDatum vd1 = new VariantDatum();
        vtsm.annotateDatum(fc, vcIndel, vd1, false);
        Assert.assertEquals(2.0, vd1.prior, 0.0, "prior");  //vd does not match -> prior does get adjusted but to default value, not the one from the feature source
        Assert.assertFalse(vd1.isKnown, "isKnown");
        Assert.assertFalse(vd1.atTruthSite, "atTruthSite");
        Assert.assertFalse(vd1.atTrainingSite, "atTrainingSite");
        Assert.assertEquals(0, vd1.consensusCount, "consensusCount");
    }

    @Test
    public void testIndelInTruthData(){
        VQSRTrainingSetManager vtsm = new VQSRTrainingSetManager();
        final String args = "--resource indels,training=true,truth=true,prior=15.0:" + TRUE_INDELS;
        TestCLP clp = new TestCLP();
        makeAndAddTrainingSetsForCLP(clp, vtsm, args);
        FeatureManager fm = new FeatureManager(clp);

        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(50).make();
        final String source= "fred";

        int start = 923978;//1	923978	122490	A	AG

        final Allele A = Allele.create("A", true);
        final Allele T = Allele.create("T");
        final List<Allele> AT = Arrays.asList(A, T);
        final double[] hetPL = MathUtils.normalizeFromRealSpace(new double[]{0.09, 0.9, 0.01});
        final Genotype atGT = new GenotypeBuilder(base).alleles(AT).AD(new int[]{10,2}).PL(hetPL).GQ(8).make();

        final List<Allele> alleles = AT;
        int stop = start + alleles.get(0).length() - 1; // alleles.contains(ATC) ? start + 3 : start;
        final Genotype genotypes = atGT;
        final String filters = null;
        final VariantContext vcSNP = new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(genotypes).filters(filters).make();

        final FeatureContext fc = new FeatureContext(fm, new SimpleInterval("1", start,	stop));
        VariantDatum vd = new VariantDatum();
        vtsm.annotateDatum(fc, vcSNP, vd, false);
        Assert.assertEquals(2.0, vd.prior, 0.0, "prior");   //vd does not match -> prior does get adjusted but to default value, not the one from the feature source
        Assert.assertFalse(vd.isKnown, "isKnown");
        Assert.assertFalse(vd.atTruthSite, "atTruthSite");
        Assert.assertFalse(vd.atTrainingSite, "atTrainingSite");
        Assert.assertEquals(0, vd.consensusCount, "consensusCount");

        final Allele AG = Allele.create("AG");
        final List<Allele> A_AG = Arrays.asList(A, AG);
        final Genotype genotypesIndel = new GenotypeBuilder(base).alleles(A_AG).AD(new int[]{10,2}).PL(hetPL).GQ(8).make();;
        final List<Allele> allelesIndel = A_AG;
        final VariantContext vcIndel = new VariantContextBuilder(source, "1", start, stop, allelesIndel).genotypes(genotypesIndel).filters(filters).make();
        VariantDatum vd1 = new VariantDatum();
        vtsm.annotateDatum(fc, vcIndel, vd1, false);
        Assert.assertEquals(15.0, vd1.prior, 0.0, "prior");
        Assert.assertFalse(vd1.isKnown, "isKnown");
        Assert.assertTrue(vd1.atTruthSite, "atTruthSite");
        Assert.assertTrue(vd1.atTrainingSite, "atTrainingSite");
        Assert.assertEquals(0, vd1.consensusCount, "consensusCount");
    }

}
