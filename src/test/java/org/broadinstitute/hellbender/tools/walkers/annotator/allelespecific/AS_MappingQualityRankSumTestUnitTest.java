package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.utils.test.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by emeryj on 8/11/17.
 */
public class AS_MappingQualityRankSumTestUnitTest extends ReducibleAnnotationBaseTest {
    private  static final Allele REF = Allele.create("T", true);
    private static final Allele ALT = Allele.create("A", false);

    private static final String sample1 = "NA1";
    private static final String sample2 = "NA2";

    private VariantContext makeVC( final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    private GATKRead makeRead(final int mq) {
        Cigar cigar = TextCigarCodec.decode("10M");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mq);
        return read;
    }

    @Override
    protected List<Annotation> getAnnotationsToUse() {
        return Collections.singletonList(new AS_MappingQualityRankSumTest());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY;
    }

    @Override
    protected String getKey() {
        return GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY;
    }

    @Test  //To be consistent with annotate raw data now produces the the ranksum, making this test incorrect
    public void testAS_MQRaw(){
        final AS_RankSumTest ann = new AS_MappingQualityRankSumTest();
        final String key1 = GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY;
        final String key2 = GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY;

        final int[] altMappingQualities = {10, 20};
        final int[] refMappingQualities = {100, 110};
        final List<GATKRead> refReads = Arrays.stream(refMappingQualities).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final List<GATKRead> altReads = Arrays.stream(altMappingQualities).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);

        final ReferenceContext ref= null;
        final VariantContext vc= makeVC(REF, ALT);

        final Map<String, Object> annotateRaw = ann.annotateRawData(ref, vc, likelihoods);
        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        MannWhitneyU.Result expectedAlt = mannWhitneyU.test(new double[]{10.0, 20.0},new double[]{100.0, 110.0}, MannWhitneyU.TestType.FIRST_DOMINATES);
        String expected = "|"+String.format("%.1f",Math.round(Math.floor((expectedAlt.getZ() )/0.1))*0.1)+",1";
        String expectedAnnotate = String.format("%.3f",expectedAlt.getZ());

        Assert.assertEquals(annotate.get(key2), expectedAnnotate);
        Assert.assertEquals(annotateRaw.get(key1), expected);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key2);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key2);
    }
}