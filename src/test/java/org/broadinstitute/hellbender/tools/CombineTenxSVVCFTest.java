package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.CombineTenxSVVCF.BreakendAdjacency;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import static org.testng.Assert.*;

public final class CombineTenxSVVCFTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testParseBreakendAllele() {
        BreakendAdjacency breakendAdjacency = CombineTenxSVVCF.parseBreakendAllele("G]17:198982]");
        Assert.assertEquals(breakendAdjacency.before, false);
        Assert.assertEquals(breakendAdjacency.revComp, true);
        Assert.assertEquals(breakendAdjacency.contig, "17");
        Assert.assertEquals(breakendAdjacency.position, 198982);
        Assert.assertEquals(breakendAdjacency.localBases, "G");
    }

    @Test(groups = "sv")
    public void testGetMaxClique() {
        final List<Allele> alleles1 = Arrays.asList(Allele.create("G", true), Allele.create("G]17:198982]"));
        final VariantContext vc1 = new VariantContextBuilder("CombineTenxSVVCFTest", "17", 150000, 150000, alleles1).make();
        final LinkedList<VariantContext> vcs1 = new LinkedList<>();
        vcs1.add(vc1);

        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("17", 1000000)));

        final Set<VariantContext> maxClique1 = CombineTenxSVVCF.getMaxClique(vcs1, 500, dictionary);
        Assert.assertEquals( maxClique1.size(), 1);
        Assert.assertTrue(maxClique1.contains(vc1));

        final VariantContext vc2 = new VariantContextBuilder("CombineTenxSVVCFTest", "17", 150900, 150900, alleles1).make();
        final LinkedList<VariantContext> vcs2 = new LinkedList<>();
        vcs2.add(vc1);
        vcs2.add(vc2);

        final Set<VariantContext> maxClique2 = CombineTenxSVVCF.getMaxClique(vcs2, 500, dictionary);
        Assert.assertEquals(maxClique2.size(), 2);
        Assert.assertTrue(maxClique2.contains(vc1));
        Assert.assertTrue(maxClique2.contains(vc2));

        final VariantContext vc3 = new VariantContextBuilder("CombineTenxSVVCFTest", "17", 151700, 151700, alleles1).make();
        final VariantContext vc4 = new VariantContextBuilder("CombineTenxSVVCFTest", "17", 151800, 151800, alleles1).make();
        final LinkedList<VariantContext> vcs3 = new LinkedList<>();
        vcs3.add(vc1);
        vcs3.add(vc2);
        vcs3.add(vc3);
        vcs3.add(vc4);

        final Set<VariantContext> maxClique3 = CombineTenxSVVCF.getMaxClique(vcs3, 500, dictionary);
        Assert.assertEquals(maxClique3.size(), 3);
        Assert.assertFalse(maxClique3.contains(vc1));
        Assert.assertTrue(maxClique3.contains(vc2));
        Assert.assertTrue(maxClique3.contains(vc3));
        Assert.assertTrue(maxClique3.contains(vc4));

        final List<Allele> alleles2 = Arrays.asList(Allele.create("G", true), Allele.create("G]17:199283]"));
        final VariantContext vc5 = new VariantContextBuilder("CombineTenxSVVCFTest", "17", 150000, 150000, alleles2).make();
        final List<Allele> alleles3 = Arrays.asList(Allele.create("C", true), Allele.create("C]17:210282]"));
        final VariantContext vc6 = new VariantContextBuilder("CombineTenxSVVCFTest", "17", 150000, 150000, alleles3).make();
        final List<Allele> alleles4 = Arrays.asList(Allele.create("A", true), Allele.create("A]17:210283]"));
        final VariantContext vc7 = new VariantContextBuilder("CombineTenxSVVCFTest", "17", 150000, 150000, alleles4).make();
        final List<Allele> alleles5 = Arrays.asList(Allele.create("A", true), Allele.create("A]17:210284]"));
        final VariantContext vc8 = new VariantContextBuilder("CombineTenxSVVCFTest", "17", 150000, 150000, alleles5).make();

        final LinkedList<VariantContext> vcs4 = new LinkedList<>();
        vcs4.add(vc1);
        vcs4.add(vc5);
        vcs4.add(vc6);
        vcs4.add(vc7);
        vcs4.add(vc8);

        final Set<VariantContext> maxClique4 = CombineTenxSVVCF.getMaxClique(vcs4, 500, dictionary);
        Assert.assertEquals(maxClique4.size(), 3);
        Assert.assertTrue(maxClique4.contains(vc6));
        Assert.assertTrue(maxClique4.contains(vc7));
        Assert.assertTrue(maxClique4.contains(vc8));

    }
}