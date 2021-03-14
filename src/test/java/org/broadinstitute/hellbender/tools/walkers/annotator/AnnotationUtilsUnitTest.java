package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.testng.Assert.*;

public class AnnotationUtilsUnitTest extends GATKBaseTest {

    @Test
    public void testGenerateMissingDataWarning() {
        final String contig = "CONTIG";
        final int start = 12345;
        final String sampleName = "sample1";

        final List<GATKRead> reads = IntStream.range(0, 10)
                .mapToObj(n -> ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"))).collect(Collectors.toList());

        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(sampleName, reads, -100.0, Allele.REF_A, Allele.ALT_T);

        final String noCallWarning = AnnotationUtils.generateMissingDataWarning(VariantContextTestUtils.makeNonRef(contig, start),
                VariantContextTestUtils.makeG(sampleName, Allele.NO_CALL, Allele.NO_CALL),
               likelihoods);
        Assert.assertTrue(noCallWarning.contains(contig));
        Assert.assertTrue(noCallWarning.contains(Integer.toString(start)));
        Assert.assertTrue(noCallWarning.contains("not called"));

        final String nullLikelihoodsWarning = AnnotationUtils.generateMissingDataWarning(VariantContextTestUtils.makeNonRef(contig, start),
                VariantContextTestUtils.makeG(sampleName, Allele.REF_A, Allele.REF_T),
                null);
        Assert.assertTrue(nullLikelihoodsWarning.contains(contig));
        Assert.assertTrue(nullLikelihoodsWarning.contains(Integer.toString(start)));
        Assert.assertTrue(nullLikelihoodsWarning.contains("alleleLikelihoodMap is null"));

        final String bothBadWarning =  AnnotationUtils.generateMissingDataWarning(VariantContextTestUtils.makeNonRef(contig, start),
                VariantContextTestUtils.makeG(sampleName, Allele.NO_CALL, Allele.NO_CALL),
                null);
        Assert.assertTrue(bothBadWarning.contains(contig));
        Assert.assertTrue(bothBadWarning.contains(Integer.toString(start)));
        Assert.assertTrue(bothBadWarning.contains("not called"));
        Assert.assertTrue(bothBadWarning.contains("alleleLikelihoodMap is null"));
    }
}