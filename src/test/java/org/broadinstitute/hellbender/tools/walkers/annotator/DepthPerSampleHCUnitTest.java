package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DepthPerSampleHCUnitTest extends GATKBaseTest {

    @Test
    public void testNoReads() {

        final Allele Aref = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AC = Arrays.asList(Aref, C);

        final ReferenceContext rc = new ReferenceContext();
        final GenotypeBuilder gb =  new GenotypeBuilder("sample", AC).DP(10).AD(new int[]{5,5});
        final Genotype g = gb.make();
        final List<GATKRead> reads = new ArrayList<>();
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods("sample", reads, -100.0, Aref, C);
        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).genotypes(Arrays.asList(g)).make();

        new DepthPerSampleHC().annotate(rc, vc, g, gb, likelihoods);
        Assert.assertEquals(gb.make().getDP(), 0);
    }
}
