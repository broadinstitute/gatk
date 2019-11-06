package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import static org.testng.Assert.*;

public class AlleleFractionGenotypeAnnotationUnitTest extends GATKBaseTest {

    @DataProvider(name = "testUsingADDataProvider")
    public Object[][] testUsingADDataProvider() {
        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AA = Arrays.asList(A, A);
        final List<Allele> AC = Arrays.asList(A, C);

        final Genotype gAC = new GenotypeBuilder("1", AC).AD(new int[]{5,5}).make();
        final Genotype gAA = new GenotypeBuilder("2", AA).AD(new int[]{10,0}).make();

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).genotypes(Arrays.asList(gAA, gAC)).make();

        return new Object[][] {
                {vc, gAC, new double[]{0.5}},
                {vc, gAA, new double[]{0.0}}
            };
    }

    @Test(dataProvider = "testUsingADDataProvider")
    public void testUsingAD(final VariantContext vc, final Genotype g, final double[] expectedAF) {
        final GenotypeBuilder gb = new GenotypeBuilder();

        new AlleleFraction().annotate(null, vc, g, gb, null);


        double[] af = (double[])gb.make().getAnyAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);
        Assert.assertEquals(af, expectedAF);
    }

    @DataProvider(name = "testUsingReadsDataProvider")
    public Object[][] testUsingReadsDataProvider() {
        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AA = Arrays.asList(A, A);
        final List<Allele> AC = Arrays.asList(A, C);

        final Genotype gAC = new GenotypeBuilder("1", AC).make();
        final Genotype gAA = new GenotypeBuilder("2", AA).make();

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).genotypes(Arrays.asList(gAA, gAC)).make();

        final byte[] qual = new byte[]{30,30,30,30};
        final String refBases = "AAAA";
        final String altBases = "AACA";

        final GATKRead refRead = ArtificialReadUtils.createArtificialRead(refBases.getBytes(), qual, "4M");
        final GATKRead altRead = ArtificialReadUtils.createArtificialRead(altBases.getBytes(), qual, "4M");

        final int nRef=7;
        final int nAlt=3;

        final List<GATKRead> sample1Reads = new ArrayList<>(Collections.nCopies(nRef, refRead));
        sample1Reads.addAll(Collections.nCopies(nAlt, altRead));

        final List<GATKRead> sample2Reads = Collections.nCopies(nRef + nAlt, refRead);

        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("1", sample1Reads, "2", sample2Reads);
        final IndexedSampleList sampleList = new IndexedSampleList(Arrays.asList("1", "2"));
        final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods = new AlleleLikelihoods<GATKRead, Allele>(sampleList, new IndexedAlleleList<>(AC), readsBySample);

        //setup likelihood matrix
        final LikelihoodMatrix<GATKRead, Allele> matrix1 = alleleLikelihoods.sampleMatrix(0);
        final LikelihoodMatrix<GATKRead, Allele> matrix2 = alleleLikelihoods.sampleMatrix(1);
        final double MATCH_LIKELIHOOD = -1.0;
        final double MISMATCH_LIKELIHOOD = -100.0;
        for (int i=0; i<nRef+nAlt; i++) {
            matrix1.set(0, i,i<nRef? MATCH_LIKELIHOOD : MISMATCH_LIKELIHOOD);
            matrix1.set(1, i, i<nRef? MISMATCH_LIKELIHOOD : MATCH_LIKELIHOOD);

            matrix2.set(0, i, MATCH_LIKELIHOOD);
            matrix2.set(1, i, MISMATCH_LIKELIHOOD);
        }

        return new Object[][] {
                {vc, gAC, alleleLikelihoods, new double[]{(double)nAlt/(double)(nAlt + nRef)}},
                {vc, gAA, alleleLikelihoods, new double[]{0.0}},
        };
    }

    @Test(dataProvider = "testUsingReadsDataProvider")
    public void testUsingRead(final VariantContext vc, final Genotype g, final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, final double[] expectedAF) {
        final GenotypeBuilder gb = new GenotypeBuilder();

        new AlleleFraction().annotate(null, vc, g, gb, alleleLikelihoods);


        double[] af = (double[])gb.make().getAnyAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);
        Assert.assertEquals(af, expectedAF);
    }
}