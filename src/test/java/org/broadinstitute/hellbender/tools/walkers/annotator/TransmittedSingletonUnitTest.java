package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.Sex;
import org.broadinstitute.hellbender.utils.samples.Trio;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class TransmittedSingletonUnitTest {

    private static final int GQ30 = 30;
    public static final String MOM = "mom";
    public static final String DAD = "dad";
    public static final String CHILD = "child";

    @DataProvider(name = "transmittedSingletons")
    public Object[][] transmittedSingletons() {
        final List<Object[]> tests = new ArrayList<>();

        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Allele noCallAllele = Allele.NO_CALL;

        final List<Allele> homRef = Arrays.asList(refAllele, refAllele);
        final List<Allele> homVar = Arrays.asList(altAllele, altAllele);
        final List<Allele> het = Arrays.asList(refAllele, altAllele);
        final List<Allele> notCall = Arrays.asList(noCallAllele, noCallAllele);

        final Genotype g00gq30 = new GenotypeBuilder("2", homRef).DP(30).PL(new int[]{0,100,100}).AD(new int[]{10, 0}).GQ(GQ30).make();
        final Genotype g01gq30 = new GenotypeBuilder("1", het).DP(30).PL(new int[]{100, 0, 100}).AD(new int[]{5, 5}).GQ(GQ30).make();
        final Genotype g11gq30 = new GenotypeBuilder("3", homVar).DP(30).PL(new int[]{100, 100, 0}).AD(new int[]{0, 10}).GQ(GQ30).make();
        final Genotype gNo = new GenotypeBuilder("4", notCall).make();

        tests.add(new Object[]{g00gq30, g00gq30, g00gq30, false, false, null});
        tests.add(new Object[]{g00gq30, g00gq30, g01gq30, false, false, null});
        tests.add(new Object[]{g00gq30, g00gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g00gq30, g00gq30, gNo, false, false, null});

        tests.add(new Object[]{g00gq30, g01gq30, g00gq30, false, true, DAD});
        tests.add(new Object[]{g00gq30, g01gq30, g01gq30, true, false, DAD});
        tests.add(new Object[]{g00gq30, g01gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g00gq30, g01gq30, gNo, false, false, null});

        tests.add(new Object[]{g00gq30, g11gq30, g00gq30, false, false, null});
        tests.add(new Object[]{g00gq30, g11gq30, g01gq30, false, false, null});
        tests.add(new Object[]{g00gq30, g11gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g00gq30, g11gq30, gNo, false, false, null});

        tests.add(new Object[]{g00gq30, gNo, g00gq30, false, false, null});
        tests.add(new Object[]{g00gq30, gNo, g01gq30, false, false, null});
        tests.add(new Object[]{g00gq30, gNo, g11gq30, false, false, null});
        tests.add(new Object[]{g00gq30, gNo, gNo, false, false, null});

        tests.add(new Object[]{g01gq30, g00gq30, g00gq30, false, true, MOM});
        tests.add(new Object[]{g01gq30, g00gq30, g01gq30, true, false, MOM});
        tests.add(new Object[]{g01gq30, g00gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g01gq30, g00gq30, gNo, false, false, null});

        tests.add(new Object[]{g01gq30, g01gq30, g00gq30, false, false, null});
        tests.add(new Object[]{g01gq30, g01gq30, g01gq30, false, false, null});
        tests.add(new Object[]{g01gq30, g01gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g01gq30, g01gq30, gNo, false, false, null});

        tests.add(new Object[]{g01gq30, g11gq30, g00gq30, false, false, null});
        tests.add(new Object[]{g01gq30, g11gq30, g01gq30, false, false, null});
        tests.add(new Object[]{g01gq30, g11gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g01gq30, g11gq30, gNo, false, false, null});

        tests.add(new Object[]{g01gq30, gNo, g00gq30, false, false, null});
        tests.add(new Object[]{g01gq30, gNo, g01gq30, false, false, null});
        tests.add(new Object[]{g01gq30, gNo, g11gq30, false, false, null});
        tests.add(new Object[]{g01gq30, gNo, gNo, false, false, null});


        tests.add(new Object[]{g11gq30, g00gq30, g00gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g00gq30, g01gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g00gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g00gq30, gNo, false, false, null});

        tests.add(new Object[]{g11gq30, g01gq30, g00gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g01gq30, g01gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g01gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g01gq30, gNo, false, false, null});

        tests.add(new Object[]{g11gq30, g11gq30, g00gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g11gq30, g01gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g11gq30, g11gq30, false, false, null});
        tests.add(new Object[]{g11gq30, g11gq30, gNo, false, false, null});

        tests.add(new Object[]{g11gq30, gNo, g00gq30, false, false, null});
        tests.add(new Object[]{g11gq30, gNo, g01gq30, false, false, null});
        tests.add(new Object[]{g11gq30, gNo, g11gq30, false, false, null});
        tests.add(new Object[]{g11gq30, gNo, gNo, false, false, null});


        tests.add(new Object[]{gNo, g00gq30, g00gq30, false, false, null});
        tests.add(new Object[]{gNo, g00gq30, g01gq30, false, false, null});
        tests.add(new Object[]{gNo, g00gq30, g11gq30, false, false, null});
        tests.add(new Object[]{gNo, g00gq30, gNo, false, false, null});

        tests.add(new Object[]{gNo, g01gq30, g00gq30, false, false, null});
        tests.add(new Object[]{gNo, g01gq30, g01gq30, false, false, null});
        tests.add(new Object[]{gNo, g01gq30, g11gq30, false, false, null});
        tests.add(new Object[]{gNo, g01gq30, gNo, false, false, null});

        tests.add(new Object[]{gNo, g11gq30, g00gq30, false, false, null});
        tests.add(new Object[]{gNo, g11gq30, g01gq30, false, false, null});
        tests.add(new Object[]{gNo, g11gq30, g11gq30, false, false, null});
        tests.add(new Object[]{gNo, g11gq30, gNo, false, false, null});

        tests.add(new Object[]{gNo, gNo, g00gq30, false, false, null});
        tests.add(new Object[]{gNo, gNo, g01gq30, false, false, null});
        tests.add(new Object[]{gNo, gNo, g11gq30, false, false, null});
        tests.add(new Object[]{gNo, gNo, gNo, false, false, null});

        final Genotype g00NoPl = new GenotypeBuilder("2", homRef).DP(30).GQ(GQ30).make();
        final Genotype g01NoPl = new GenotypeBuilder("1", het).DP(30).GQ(GQ30).make();
        final Genotype g11NoPl = new GenotypeBuilder("3", homVar).DP(30).GQ(GQ30).make();

        tests.add(new Object[]{g00NoPl, g01NoPl, g11NoPl, false, false, null});
        tests.add(new Object[]{g00NoPl, g00NoPl, g01NoPl, false, false, null});
        tests.add(new Object[]{g00NoPl, g01NoPl, g01NoPl, true, false, DAD});
        tests.add(new Object[]{g00NoPl, g01NoPl, g00NoPl, false, true, DAD});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "transmittedSingletons")
    public void testUsingVC(final Genotype gMom, final Genotype gDad, final Genotype gChild, final boolean expectedTransmittedSingleton,
                            final boolean expectedNonTransmittedSingleton, final String matchingParentSample) {

        final Sample sMom = new Sample(MOM, "fam", null, null, Sex.FEMALE);
        final Sample sDad = new Sample(DAD,     sMom.getFamilyID(), null, null, Sex.MALE);
        final Sample sChild = new Sample(CHILD, sMom.getFamilyID(), sDad.getID(), sMom.getID(), Sex.MALE);

        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Allele noCallAllele = Allele.NO_CALL;

        final ArrayList<Genotype> genotypes = new ArrayList<>(Arrays.<Genotype>asList(gMom, gDad, gChild));
        final Map<String, Integer> sampleNameToOffset= new LinkedHashMap<>(3);
        sampleNameToOffset.put(MOM, 0);
        sampleNameToOffset.put(DAD, 1);
        sampleNameToOffset.put(CHILD, 2);
        final List<String> sampleNamesInOrder=Arrays.asList(MOM, DAD, CHILD);
        final GenotypesContext context = GenotypesContext.create(genotypes, sampleNameToOffset, sampleNamesInOrder);

        final Collection<Allele> alleles= new LinkedHashSet<>(2);
        alleles.addAll(gMom.getAlleles());
        alleles.addAll(gDad.getAlleles());
        alleles.addAll(gChild.getAlleles());
        alleles.add(refAllele);
        alleles.add(altAllele);
        alleles.remove(noCallAllele);
        final int ac = (int) context.stream().mapToLong(g-> g.getAlleles().stream().filter(a->a.isCalled()&&a.isNonReference()).count()).sum();
        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(context).make();
        vc.getCommonInfo().putAttribute(VCFConstants.ALLELE_COUNT_KEY, ac);

        final Trio trio = new Trio(sMom, sDad, sChild);
        final InfoFieldAnnotation ann = new TransmittedSingleton(Collections.singleton(trio));
        final Map<String, Object> result = ann.annotate(null, vc, null);
        if (expectedTransmittedSingleton) {
            Assert.assertEquals(result.get(GATKVCFConstants.TRANSMITTED_SINGLETON), Collections.singletonList(matchingParentSample));
        } else if (expectedNonTransmittedSingleton) {
            Assert.assertEquals(result.get(GATKVCFConstants.NON_TRANSMITTED_SINGLETON), Collections.singletonList(matchingParentSample));
        } else {
            Assert.assertTrue(result.isEmpty());
        }

        Assert.assertEquals(ann.getDescriptions(), Arrays.asList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.TRANSMITTED_SINGLETON),
                GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NON_TRANSMITTED_SINGLETON)));
        Assert.assertEquals(ann.getKeyNames(), Arrays.asList(GATKVCFConstants.TRANSMITTED_SINGLETON, GATKVCFConstants.NON_TRANSMITTED_SINGLETON));

        final Map<String, Object> resultNoTrio = new TransmittedSingleton(Collections.emptySet()).annotate(null, vc, null);
        Assert.assertTrue(resultNoTrio.isEmpty());
    }

    // tests quad family which won't work for transmitted singletons
    @Test
    public void testNonTrioFamilies() {
        final Sample sMom = new Sample(MOM, "fam", null, null, Sex.FEMALE);
        final Sample sDad = new Sample(DAD,     sMom.getFamilyID(), null, null, Sex.MALE);
        final Sample sSib1 = new Sample(CHILD, sMom.getFamilyID(), sDad.getID(), sMom.getID(), Sex.MALE);
        final Sample sSib2 = new Sample(CHILD + "2", sMom.getFamilyID(), sDad.getID(), sMom.getID(), Sex.FEMALE);

        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final List<Allele> homRef = Arrays.asList(refAllele, refAllele);
        final List<Allele> het = Arrays.asList(refAllele, altAllele);

        final Genotype gMom = new GenotypeBuilder("2", homRef).DP(30).PL(new int[]{0,100,100}).AD(new int[]{10, 0}).GQ(GQ30).make();
        final Genotype gDad = new GenotypeBuilder("1", het).DP(30).PL(new int[]{100, 0, 100}).AD(new int[]{5, 5}).GQ(GQ30).make();
        final Genotype gSib1 = new GenotypeBuilder("3", het).DP(30).PL(new int[]{100, 100, 0}).AD(new int[]{0, 10}).GQ(GQ30).make();
        final Genotype gSib2 = new GenotypeBuilder("4", het).DP(30).PL(new int[]{100, 0, 100}).AD(new int[]{5, 5}).GQ(GQ30).make();

        final ArrayList<Genotype> genotypes = new ArrayList<>(Arrays.<Genotype>asList(gMom, gDad, gSib1, gSib2));
        final Map<String, Integer> sampleNameToOffset= new LinkedHashMap<>(3);
        sampleNameToOffset.put(MOM, 0);
        sampleNameToOffset.put(DAD, 1);
        sampleNameToOffset.put(CHILD, 2);
        sampleNameToOffset.put(CHILD + "2", 3);
        final List<String> sampleNamesInOrder=Arrays.asList(MOM, DAD, CHILD, CHILD+"2");
        final GenotypesContext context = GenotypesContext.create(genotypes, sampleNameToOffset, sampleNamesInOrder);

        final Collection<Allele> alleles= new LinkedHashSet<>(2);
        alleles.addAll(gMom.getAlleles());
        alleles.addAll(gDad.getAlleles());
        alleles.addAll(gSib1.getAlleles());
        alleles.addAll(gSib2.getAlleles());
        alleles.add(refAllele);
        alleles.add(altAllele);
        final int ac = (int) context.stream().mapToLong(g-> g.getAlleles().stream().filter(a->a.isCalled()&&a.isNonReference()).count()).sum();
        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(context).make();
        vc.getCommonInfo().putAttribute(VCFConstants.ALLELE_COUNT_KEY, ac);

        final Trio trio1 = new Trio(sMom, sDad, sSib1);
        final Trio trio2 = new Trio(sMom, sDad, sSib2);

        final InfoFieldAnnotation ann = new TransmittedSingleton(new HashSet<>(Arrays.asList(trio1, trio2)));
        final Map<String, Object> result = ann.annotate(null, vc, null);

        //quads are not implemented yet so we skip over them
        Assert.assertEquals(result.size(), 0);
    }

}
