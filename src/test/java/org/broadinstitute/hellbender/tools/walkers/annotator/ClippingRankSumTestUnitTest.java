package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class ClippingRankSumTestUnitTest {
    private static final Allele REF = Allele.create("T", true);
    private static final Allele ALT = Allele.create("A", false);

    private static final String SAMPLE_1 = "NA1";
    private static final String SAMPLE_2 = "NA2";

    private VariantContext makeVC( final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // SAMPLE_1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // SAMPLE_2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(SAMPLE_2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    private static GATKRead makeRead(final int hardClip) {
        Cigar cigar = hardClip == 0 ? TextCigarCodec.decode("10M") : TextCigarCodec.decode("10M" + hardClip + "H");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(30);
        return read;
    }

    @Test
    public void testClipping(){
        final int[] refHardClips = {10, 0};
        final int[] altHardClips = {1, 2};
        final List<GATKRead> refReads = Arrays.stream(refHardClips).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final List<GATKRead> altReads = Arrays.stream(altHardClips).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE_1, refReads, altReads, -100.0, -100.0, REF, ALT);

        final ReferenceContext ref= null;
        final VariantContext vc= makeVC(REF, ALT);
        final InfoFieldAnnotation ann = new ClippingRankSumTest();
        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);

        final double zScore = mannWhitneyU.test(new double[]{altHardClips[0], altHardClips[1]}, new double[]{refHardClips[0], refHardClips[1]}, MannWhitneyU.TestType.FIRST_DOMINATES).getZ();
        final String zScoreStr = String.format("%.3f", zScore);
        Assert.assertEquals(annotate.get(GATKVCFConstants.CLIPPING_RANK_SUM_KEY), zScoreStr);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), GATKVCFConstants.CLIPPING_RANK_SUM_KEY);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), GATKVCFConstants.CLIPPING_RANK_SUM_KEY);


    }
}
