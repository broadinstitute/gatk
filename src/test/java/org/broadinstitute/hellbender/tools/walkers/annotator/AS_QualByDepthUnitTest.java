package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class AS_QualByDepthUnitTest extends BaseTest {

    @Test
    public void testAnnotate() throws Exception {
        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");
        final Allele G = Allele.create("G");

        final List<Allele> AC = Arrays.asList(A, C);
        final int readDepth = 20;
        final String sample1 = "sample1";
        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(sample1, AC).DP(dpDepth).make();

        final double log10PError = -5;

        final int n1A= readDepth;
        for (int i = 0; i < n1A; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "n1A_" + i);
            read.setMappingQuality(20);
            map.add(read, A, -1.0);
            map.add(read, C, -100.0);  //try to fool it - add another likelihood to same read
            map.add(read, G, -1000.0);  //and a third one
        }

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = Collections.singletonMap(sample1, map);
        final Map<String, Object> annotatedMap = new AS_QualByDepth().annotate(null, vc, perReadAlleleLikelihoodMap);
        Assert.assertNull(annotatedMap);

        final Map<String, Object> annotatedMapRaw = new AS_QualByDepth().annotateRawData(null, vc, perReadAlleleLikelihoodMap);
        Assert.assertNull(annotatedMapRaw);

    }

    @Test
    public void testDescriptions() throws Exception {
        final AS_QualByDepth cov = new AS_QualByDepth();
        Assert.assertEquals(cov.getRawKeyName(), GATKVCFConstants.AS_QUAL_KEY);
        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY);
    }

}
