package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;


public final class SampleListUnitTest extends GATKBaseTest {

    private static final int GQ30 = 30;

    @Test
    public void testUsingVC() {

        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Allele noCallAllele = Allele.NO_CALL;

        final List<Allele> homRef = Arrays.asList(refAllele, refAllele);
        final List<Allele> homVar = Arrays.asList(altAllele, altAllele);
        final List<Allele> het = Arrays.asList(refAllele, altAllele);
        final List<Allele> notCall = Arrays.asList(noCallAllele, noCallAllele);

        final Genotype g00Mom = new GenotypeBuilder("mom", homRef).DP(10).PL(new int[]{0,100,100}).AD(new int[]{10, 0}).GQ(GQ30).make();
        final Genotype g01Dad = new GenotypeBuilder("dad", het).DP(10).PL(new int[]{100, 0, 100}).AD(new int[]{5, 5}).GQ(GQ30).make();
        final Genotype g01Child = new GenotypeBuilder("child", het).DP(10).PL(new int[]{100, 0, 100}).AD(new int[]{5, 5}).GQ(GQ30).make();

        final Genotype g00Dad = new GenotypeBuilder("dad", homRef).DP(10).PL(new int[]{0,100,100}).AD(new int[]{10,0}).GQ(GQ30).make();
        final Genotype g00Child = new GenotypeBuilder("child", homRef).DP(10).PL(new int[]{0,100,100}).AD(new int[]{10,0}).GQ(GQ30).make();


        final ArrayList<Genotype> genotypesHomHetHet = new ArrayList<>(Arrays.<Genotype>asList(g01Child, g01Dad, g00Mom));
        final Map<String, Integer> sampleNameToOffset= new LinkedHashMap<>(3);
        sampleNameToOffset.put("child", 0);
        sampleNameToOffset.put("dad", 1);
        sampleNameToOffset.put("mom", 2);
        final List<String> sampleNamesInOrder=Arrays.asList("child", "dad", "mom");
        final GenotypesContext contextHomHetHet = GenotypesContext.create(genotypesHomHetHet, sampleNameToOffset, sampleNamesInOrder);

        final Collection<Allele> alleles= new LinkedHashSet<>(2);
        alleles.addAll(g00Mom.getAlleles());
        alleles.addAll(g01Dad.getAlleles());
        alleles.addAll(g01Child.getAlleles());
        alleles.add(refAllele);
        alleles.add(altAllele);
        alleles.remove(noCallAllele);

        final VariantContext vcHomHetHet = new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(contextHomHetHet).make();

        final InfoFieldAnnotation ann = new SampleList();
        final Map<String, Object> result = ann.annotate(null, vcHomHetHet, null);
        Assert.assertEquals(result.get(GATKVCFConstants.SAMPLE_LIST_KEY), "child,dad");

        final VariantContext vcNoCalls = new VariantContextBuilder("test", "20", 10, 10, alleles).make();
        final Map<String, Object> resultNoCalls = ann.annotate(null, vcNoCalls, null);
        Assert.assertTrue(resultNoCalls.isEmpty());

        final ArrayList<Genotype> genotypesHomHomHom = new ArrayList<>(Arrays.<Genotype>asList(g00Child, g00Dad, g00Mom));
        final GenotypesContext contextHomHomHom = GenotypesContext.create(genotypesHomHomHom, sampleNameToOffset, sampleNamesInOrder);

        final VariantContext vcHomHomHom = new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(contextHomHomHom).make();
        final Map<String, Object> resultHomHomHom = ann.annotate(null, vcHomHomHom, null);
        Assert.assertTrue(resultHomHomHom.isEmpty());


        Assert.assertEquals(ann.getDescriptions(), Arrays.asList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SAMPLE_LIST_KEY)));
        Assert.assertEquals(ann.getKeyNames(), Arrays.asList(GATKVCFConstants.SAMPLE_LIST_KEY));
    }

    @Test
    public void testEmptyIfNoGenotypes() throws Exception {
        final SampleList ann = new SampleList();
        final Map<String, Object> annotate = ann.annotate(null, when(mock(VariantContext.class).getGenotypesOrderedByName()).thenReturn(Collections.<Genotype>emptyList()).getMock(), null);
        Assert.assertTrue(annotate.isEmpty());
    }
}
