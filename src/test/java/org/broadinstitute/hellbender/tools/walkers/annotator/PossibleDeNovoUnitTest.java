package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.Sex;
import org.broadinstitute.hellbender.utils.samples.Trio;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class PossibleDeNovoUnitTest extends GATKBaseTest {

    private static final int GQ30 = 30;
    private static final int GQ45 = 45;

    @DataProvider(name = "deNovo")
    public Object[][] deNovo() {
        final List<Object[]> tests = new ArrayList<>();

        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Allele noCallAllele = Allele.NO_CALL;

        final List<Allele> homRef = Arrays.asList(refAllele, refAllele);
        final List<Allele> homVar = Arrays.asList(altAllele, altAllele);
        final List<Allele> het = Arrays.asList(refAllele, altAllele);
        final List<Allele> notCall = Arrays.asList(noCallAllele, noCallAllele);

        final Genotype g00gq30 = new GenotypeBuilder("2", homRef).DP(10).PL(new int[]{0,100,100}).AD(new int[]{10, 0}).GQ(GQ30).make();
        final Genotype g01gq30 = new GenotypeBuilder("1", het).DP(10).PL(new int[]{100, 0, 100}).AD(new int[]{5, 5}).GQ(GQ30).make();
        final Genotype g11gq30 = new GenotypeBuilder("3", homVar).DP(10).PL(new int[]{100, 100, 0}).AD(new int[]{0, 10}).GQ(GQ30).make();
        final Genotype gNo = new GenotypeBuilder("4", notCall).make();

        tests.add(new Object[]{g00gq30, g00gq30, g00gq30, false, false});
        //test case where there is a possible DeNovo, but depth and parent GQs are lower than set thresholds (20 and 40 respectively)
        tests.add(new Object[]{g00gq30, g00gq30, g01gq30, true, false});
        tests.add(new Object[]{g00gq30, g00gq30, g11gq30, false, false});
        tests.add(new Object[]{g00gq30, g00gq30, gNo, false, false});

        tests.add(new Object[]{g00gq30, g01gq30, g00gq30, false, false});
        tests.add(new Object[]{g00gq30, g01gq30, g01gq30, false, false});
        tests.add(new Object[]{g00gq30, g01gq30, g11gq30, false, false});
        tests.add(new Object[]{g00gq30, g01gq30, gNo, false, false});

        tests.add(new Object[]{g00gq30, g11gq30, g00gq30, false, false});
        tests.add(new Object[]{g00gq30, g11gq30, g01gq30, false, false});
        tests.add(new Object[]{g00gq30, g11gq30, g11gq30, false, false});
        tests.add(new Object[]{g00gq30, g11gq30, gNo, false, false});

        tests.add(new Object[]{g00gq30, gNo, g00gq30, false, false});
        tests.add(new Object[]{g00gq30, gNo, g01gq30, false, false});
        tests.add(new Object[]{g00gq30, gNo, g11gq30, false, false});
        tests.add(new Object[]{g00gq30, gNo, gNo, false, false});

        tests.add(new Object[]{g01gq30, g00gq30, g00gq30, false, false});
        tests.add(new Object[]{g01gq30, g00gq30, g01gq30, false, false});
        tests.add(new Object[]{g01gq30, g00gq30, g11gq30, false, false});
        tests.add(new Object[]{g01gq30, g00gq30, gNo, false, false});

        tests.add(new Object[]{g01gq30, g01gq30, g00gq30, false, false});
        tests.add(new Object[]{g01gq30, g01gq30, g01gq30, false, false});
        tests.add(new Object[]{g01gq30, g01gq30, g11gq30, false, false});
        tests.add(new Object[]{g01gq30, g01gq30, gNo, false, false});

        tests.add(new Object[]{g01gq30, g11gq30, g00gq30, false, false});
        tests.add(new Object[]{g01gq30, g11gq30, g01gq30, false, false});
        tests.add(new Object[]{g01gq30, g11gq30, g11gq30, false, false});
        tests.add(new Object[]{g01gq30, g11gq30, gNo, false, false});

        tests.add(new Object[]{g01gq30, gNo, g00gq30, false, false});
        tests.add(new Object[]{g01gq30, gNo, g01gq30, false, false});
        tests.add(new Object[]{g01gq30, gNo, g11gq30, false, false});
        tests.add(new Object[]{g01gq30, gNo, gNo, false, false});


        tests.add(new Object[]{g11gq30, g00gq30, g00gq30, false, false});
        tests.add(new Object[]{g11gq30, g00gq30, g01gq30, false, false});
        tests.add(new Object[]{g11gq30, g00gq30, g11gq30, false, false});
        tests.add(new Object[]{g11gq30, g00gq30, gNo, false, false});

        tests.add(new Object[]{g11gq30, g01gq30, g00gq30, false, false});
        tests.add(new Object[]{g11gq30, g01gq30, g01gq30, false, false});
        tests.add(new Object[]{g11gq30, g01gq30, g11gq30, false, false});
        tests.add(new Object[]{g11gq30, g01gq30, gNo, false, false});

        tests.add(new Object[]{g11gq30, g11gq30, g00gq30, false, false});
        tests.add(new Object[]{g11gq30, g11gq30, g01gq30, false, false});
        tests.add(new Object[]{g11gq30, g11gq30, g11gq30, false, false});
        tests.add(new Object[]{g11gq30, g11gq30, gNo, false, false});

        tests.add(new Object[]{g11gq30, gNo, g00gq30, false, false});
        tests.add(new Object[]{g11gq30, gNo, g01gq30, false, false});
        tests.add(new Object[]{g11gq30, gNo, g11gq30, false, false});
        tests.add(new Object[]{g11gq30, gNo, gNo, false, false});


        tests.add(new Object[]{gNo, g00gq30, g00gq30, false, false});
        tests.add(new Object[]{gNo, g00gq30, g01gq30, false, false});
        tests.add(new Object[]{gNo, g00gq30, g11gq30, false, false});
        tests.add(new Object[]{gNo, g00gq30, gNo, false, false});

        tests.add(new Object[]{gNo, g01gq30, g00gq30, false, false});
        tests.add(new Object[]{gNo, g01gq30, g01gq30, false, false});
        tests.add(new Object[]{gNo, g01gq30, g11gq30, false, false});
        tests.add(new Object[]{gNo, g01gq30, gNo, false, false});

        tests.add(new Object[]{gNo, g11gq30, g00gq30, false, false});
        tests.add(new Object[]{gNo, g11gq30, g01gq30, false, false});
        tests.add(new Object[]{gNo, g11gq30, g11gq30, false, false});
        tests.add(new Object[]{gNo, g11gq30, gNo, false, false});

        tests.add(new Object[]{gNo, gNo, g00gq30, false, false});
        tests.add(new Object[]{gNo, gNo, g01gq30, false, false});
        tests.add(new Object[]{gNo, gNo, g11gq30, false, false});
        tests.add(new Object[]{gNo, gNo, gNo, false, false});

        final Genotype g00NoPl = new GenotypeBuilder("2", homRef).DP(10).GQ(GQ30).make();
        final Genotype g01NoPl = new GenotypeBuilder("1", het).DP(10).GQ(GQ30).make();
        final Genotype g11NoPl = new GenotypeBuilder("3", homVar).DP(10).GQ(GQ30).make();

        tests.add(new Object[]{g00NoPl, g01NoPl, g11NoPl, false, false});
        tests.add(new Object[]{g00NoPl, g00NoPl, g01NoPl, true, false});

        final Genotype g00gq45 = new GenotypeBuilder("2", homRef).DP(50).PL(new int[]{0,100,100}).AD(new int[]{10, 0}).GQ(GQ45).make();
        final Genotype g01gq45 = new GenotypeBuilder("1", het).DP(50).PL(new int[]{100, 0, 100}).AD(new int[]{5, 5}).GQ(GQ45).make();

        //test case when there is a possible deNovo (het child with hom ref parents) and GQ is high enough to pass parent GQ threshold of 40 and depth threshold of 20
        tests.add(new Object[]{g00gq45, g00gq45, g01gq45, true, true});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "deNovo")
    public void testUsingVC(final Genotype gMom, final Genotype gDad, final Genotype gChild, final boolean expectedDeNovo, final boolean expectedDeNovoWithExtraThresholds) {

        final Sample sMom = new Sample("mom", "fam", null, null, Sex.FEMALE);
        final Sample sDad = new Sample("dad",     sMom.getFamilyID(), null, null, Sex.MALE);
        final Sample sChild = new Sample("child", sMom.getFamilyID(), sDad.getID(), sMom.getID(), Sex.MALE);

        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Allele noCallAllele = Allele.NO_CALL;

        final ArrayList<Genotype> genotypes = new ArrayList<>(Arrays.<Genotype>asList(gMom, gDad, gChild));
        final Map<String, Integer> sampleNameToOffset= new LinkedHashMap<>(3);
        sampleNameToOffset.put("mom", 0);
        sampleNameToOffset.put("dad", 1);
        sampleNameToOffset.put("child", 2);
        final List<String> sampleNamesInOrder=Arrays.asList("mom", "dad", "child");
        final GenotypesContext context = GenotypesContext.create(genotypes, sampleNameToOffset, sampleNamesInOrder);

        final Collection<Allele> alleles= new LinkedHashSet<>(2);
        alleles.addAll(gMom.getAlleles());
        alleles.addAll(gDad.getAlleles());
        alleles.addAll(gChild.getAlleles());
        alleles.add(refAllele);
        alleles.add(altAllele);
        alleles.remove(noCallAllele);
        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(context).make();

        final Trio trio = new Trio(sMom, sDad, sChild);
        final InfoFieldAnnotation ann = new PossibleDeNovo(Collections.singleton(trio), 0);
        final InfoFieldAnnotation hiQualAnn = new PossibleDeNovo(Collections.singleton(trio), 0, 40, 20);
        final Map<String, Object> result = ann.annotate(null, vc, null);
        final Map<String, Object> hiQualResult = hiQualAnn.annotate(null, vc, null);
        if (expectedDeNovo) {
            Assert.assertEquals(result.get(GATKVCFConstants.HI_CONF_DENOVO_KEY), Collections.singletonList(sChild.getID()));
            Assert.assertFalse(result.containsKey(GATKVCFConstants.LO_CONF_DENOVO_KEY));
        }
        if (expectedDeNovoWithExtraThresholds) {
            Assert.assertEquals(hiQualResult.get(GATKVCFConstants.HI_CONF_DENOVO_KEY), Collections.singletonList(sChild.getID()));
            Assert.assertFalse(result.containsKey(GATKVCFConstants.LO_CONF_DENOVO_KEY));
        }
        if (!expectedDeNovo && !expectedDeNovoWithExtraThresholds) {
            Assert.assertTrue(result.isEmpty());
        }

        Assert.assertEquals(ann.getDescriptions(), Arrays.asList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.HI_CONF_DENOVO_KEY),
                                                                               GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.LO_CONF_DENOVO_KEY)));
        Assert.assertEquals(ann.getKeyNames(), Arrays.asList(GATKVCFConstants.HI_CONF_DENOVO_KEY, GATKVCFConstants.LO_CONF_DENOVO_KEY));

        final Map<String, Object> resultNoTrio = new PossibleDeNovo(Collections.emptySet(), 0).annotate(null, vc, null);
        Assert.assertTrue(resultNoTrio.isEmpty());
    }
}
