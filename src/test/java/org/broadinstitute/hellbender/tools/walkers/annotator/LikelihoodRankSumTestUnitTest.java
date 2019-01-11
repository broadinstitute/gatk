package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

public final class LikelihoodRankSumTestUnitTest extends GATKBaseTest {
    final Allele REF = Allele.create("T", true);
    final Allele ALT = Allele.create("A", false);

    private static final String sample1 = "NA1";
    private static final String sample2 = "NA2";
    private static final String CONTIG = "chr1";

    private VariantContext makeVC( final String contig, final long position, final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr(contig).start(position).stop(position).genotypes(testGC).make();
    }

    private GATKRead makeRead() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
        read.setMappingQuality(30);
        read.setPosition(CONTIG, 1);
        return read;
    }

    @Test
    public void testReadPos(){
        final double[] altBestAlleleLL = {-1.0, -2.0};
        final double[] refBestAlleleLL = {-5.0, -7.0};
        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        final List<GATKRead> refReads = Arrays.asList(makeRead(), makeRead());
        final List<GATKRead> altReads = Arrays.asList(makeRead(), makeRead());

        // first two reads are ref, last two are alt, "wrong" likelihoods are -100
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);

        // modify "good" likelihoods manually
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        matrix.set(0, 0, refBestAlleleLL[0]);
        matrix.set(0, 1, refBestAlleleLL[1]);
        matrix.set(1, 2, altBestAlleleLL[0]);
        matrix.set(1, 3, altBestAlleleLL[1]);
        
        final InfoFieldAnnotation ann = new LikelihoodRankSumTest();
        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), GATKVCFConstants.LIKELIHOOD_RANK_SUM_KEY);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), GATKVCFConstants.LIKELIHOOD_RANK_SUM_KEY);

        final ReferenceContext ref= null;

        final long position = 5L;  //middle of the read
        final VariantContext vc= makeVC(CONTIG, position, REF, ALT);

        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);
        final double zScore = mannWhitneyU.test(new double[]{altBestAlleleLL[0], altBestAlleleLL[1]}, new double[]{refBestAlleleLL[0], refBestAlleleLL[1]}, MannWhitneyU.TestType.FIRST_DOMINATES).getZ();
        final String zScoreStr= String.format("%.3f", zScore);
        Assert.assertEquals(annotate.get(GATKVCFConstants.LIKELIHOOD_RANK_SUM_KEY), zScoreStr);
    }

}
